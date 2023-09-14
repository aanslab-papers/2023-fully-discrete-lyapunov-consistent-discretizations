"""SSDC I/O utilities."""
# pylint: disable=invalid-name
# pylint: disable=missing-docstring
# pylint: disable=consider-using-f-string
# pylint: disable=too-many-public-methods
# pylint: disable=useless-object-inheritance
# pylint: disable=too-many-instance-attributes
from __future__ import print_function
from __future__ import absolute_import
from collections import OrderedDict
import numpy as np


_PRECISION = {
    'single': {'real': '>f4', 'complex': '>c8'},
    'double': {'real': '>f8', 'complex': '>c16'},
}

_INDEXTYPE = {
    '32bit': '>i4',
    '64bit': '>i8',
}


def _fread(stream, dtype, count=None):
    assert stream is not None
    assert dtype is not None
    dtype = np.dtype(dtype)
    size = (count if count is not None else 1)
    data = np.fromfile(stream, dtype, size)
    if data.size < size:
        raise EOFError
    data = (data if count is not None else data[0])
    return data.astype(dtype.newbyteorder('='))


def _fwrite(stream, dtype, data):
    assert stream is not None
    assert dtype is not None
    assert data is not None
    dtype = np.dtype(dtype)
    np.asarray(data, dtype).tofile(stream)


class IO(object):

    _registry = {}

    @classmethod
    def register(cls, klass):
        cls._registry[klass.CLASSID] = klass
        return klass

    def __init__(self, precision='double', scalar='real', indices='32bit'):
        self._tp_int = np.dtype(_INDEXTYPE[indices])
        self._tp_real = np.dtype(_PRECISION[precision]['real'])
        self._tp_complex = np.dtype(_PRECISION[precision]['complex'])
        self._tp_scalar = np.dtype(_PRECISION[precision][scalar])
        self.stream = None

    def __del__(self):
        self.close()

    def open(self, file, mode='r'):
        # pylint: disable=consider-using-with
        # pylint: disable=unspecified-encoding
        assert self.stream is None
        self.stream = open(file, mode+'b')
        return self

    def close(self):
        if self.stream is not None:
            self.stream.close()
            self.stream = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def _fread(self, dtype, count=None):
        return _fread(self.stream, dtype, count)

    def _fwrite(self, dtype, data):
        return _fwrite(self.stream, dtype, data)

    def load_size(self, count=None):
        size = self.load_int(count)
        assert np.min(size) >= 0
        return size

    def dump_size(self, obj):
        assert np.min(obj) >= 0
        self.dump_int(obj)

    def load_int(self, count=None):
        return self._fread(self._tp_int, count)

    def dump_int(self, obj):
        return self._fwrite(self._tp_int, obj)

    def load_real(self, count=None):
        return self._fread(self._tp_real, count)

    def dump_real(self, obj):
        return self._fwrite(self._tp_real, obj)

    def load_complex(self, count=None):
        return self._fread(self._tp_complex, count)

    def dump_complex(self, obj):
        return self._fwrite(self._tp_complex, obj)

    def load_scalar(self, count=None):
        return self._fread(self._tp_scalar, count)

    def dump_scalar(self, obj):
        return self._fwrite(self._tp_scalar, obj)

    def load_string(self, size):
        assert size >= 0
        dtype = np.dtype('S{:d}'.format(size))
        return self._fread(dtype).decode('utf-8')

    def dump_string(self, obj, size):
        assert size >= 0
        dtype = np.dtype('S{:d}'.format(size))
        self._fwrite(dtype, obj.encode('utf-8'))

    def load_cid(self, cls=None):
        cid = int(self.load_int())
        if cls is not None:
            assert cls.CLASSID is not None
            assert cls.CLASSID == cid
        return cid

    def dump_cid(self, cls):
        assert cls.CLASSID is not None
        self.dump_int(cls.CLASSID)

    def load_cls(self, cls):
        return cls.load(self)

    def dump_cls(self, cls, obj):
        cls.dump(self, obj)

    def peek(self):
        off = self.stream.tell()
        try:
            cid = self.load_cid()
            return self._registry[cid]
        finally:
            self.stream.seek(off)

    def load(self, cls=None):
        if cls is None:
            cls = self.peek()
        return self.load_cls(cls)

    def dump(self, cls, obj):
        assert cls is not None
        return self.dump_cls(cls, obj)


@IO.register
class IS(object):

    CLASSID = 1211218

    @classmethod
    def load(cls, io):
        io.load_cid(cls)
        size = io.load_size()
        data = io.load_int(size)
        return data

    @classmethod
    def dump(cls, io, obj):
        io.dump_cid(cls)
        io.dump_size(obj.size)
        io.dump_int(obj)


@IO.register
class Vec(object):

    CLASSID = 1211214

    @classmethod
    def load(cls, io):
        io.load_cid(cls)
        size = io.load_size()
        data = io.load_scalar(size)
        return data

    @classmethod
    def dump(cls, io, obj):
        io.dump_cid(cls)
        io.dump_size(obj.size)
        io.dump_scalar(obj)


@IO.register
class Mat(object):

    CLASSID = 1211216

    @classmethod
    def load(cls, io):
        io.load_cid(cls)
        M, N = io.load_size(2)
        nnz = io.load_int()
        assert nnz >= 0 or nnz == -1
        sparse = (nnz >= 0)
        if sparse:
            rnz = io.load_size(M)
            I = np.empty(M+1, rnz.dtype)
            I[0] = 0
            np.cumsum(rnz, out=I[1:])
            assert I[-1] == nnz
            J = io.load_int(nnz)
            assert np.min(J) >= 0
            assert np.max(J) < N
            V = io.load_scalar(nnz)
            data = ((M, N), (I, J, V))
        else:
            V = io.load_scalar(M*N)
            data = V.reshape(M, N)
        return data

    @classmethod
    def dump(cls, io, obj):
        try:
            (M, N), (I, J, V) = obj
            sparse = True
        except (TypeError, ValueError):
            (M, N) = obj.shape, obj
            (I, J, V) = None, None, obj
            sparse = False
        if sparse:
            nnz = J.size
            rnz = np.diff(I)
            assert I.size == M + 1
            assert I[0] == 0
            assert I[-1] == nnz
            assert np.min(J) >= 0
            assert np.max(J) < N
            assert V.size == nnz
        else:
            nnz = -1
            rnz = None
        io.dump_cid(cls)
        io.dump_size([M, N])
        io.dump_int(nnz)
        if sparse:
            io.dump_size(rnz)
            io.dump_int(J)
        io.dump_scalar(V)


@IO.register
class DMPlex(object):

    CLASSID = 1211221

    def __init__(self):
        self.dim = None
        self.point = None
        self.sizes = None
        self.cones = None
        self.ornts = None
        self.reftree = None
        self.parents = None
        self.childid = None
        self.nsd = None
        self.coords = None
        self.labels = None

    @classmethod
    def load(cls, io):
        dm = cls()
        io.load_cid(cls)
        dm.dim = io.load_int()
        assert 0 <= dm.dim <= 3

        dm.point = io.load_cls(IS)
        dm.sizes = io.load_cls(IS)
        dm.cones = io.load_cls(IS)
        dm.ornts = io.load_cls(IS)

        mask = int(io.load_int())
        if mask & 1:
            dm.reftree = io.load_cls(DMPlex)
        if mask & 2:
            keys = io.load_cls(IS)
            vals = io.load_cls(IS)
            dm.parents = (keys, vals)
        if mask & 1 and mask & 2:
            dm.childid = io.load_cls(IS)

        dm.nsd = io.load_int()
        assert -1 <= dm.nsd <= 3
        mask = int(io.load_int())
        if mask:
            dm.coords = [None] * 2
            for i in range(2):
                if mask & (1 << i):
                    dm.coords[i] = cls._load_coords(io)
            dm.coords = tuple(dm.coords)

        dm.labels = OrderedDict()
        count = io.load_size()
        for _ in range(count):
            label = cls._load_label(io)
            dm.labels[label.name] = label

        return dm

    @classmethod
    def dump(cls, io, obj):
        dm = obj
        io.dump_cid(cls)
        io.dump_int(dm.dim)

        io.dump_cls(IS, dm.point)
        io.dump_cls(IS, dm.sizes)
        io.dump_cls(IS, dm.cones)
        io.dump_cls(IS, dm.ornts)

        mask = 0
        if dm.reftree is not None:
            mask |= 1
        if dm.parents is not None:
            mask |= 2
        io.dump_int(mask)
        if mask & 1:
            io.dump_cls(DMPlex, dm.reftree)
        if mask & 2:
            io.dump_cls(IS, dm.parents[0])
            io.dump_cls(IS, dm.parents[1])
        if mask & 1 and mask & 2:
            io.dump_cls(IS, dm.childid)

        io.dump_int(dm.nsd)
        mask = 0
        if dm.coords is not None:
            for i, coords in enumerate(dm.coords):
                if coords is not None:
                    mask |= (1 << i)
        io.dump_int(mask)
        if dm.coords is not None:
            for coords in dm.coords:
                if coords is not None:
                    cls._dump_coords(io, coords)

        io.dump_size(len(dm.labels))
        for label in dm.labels.values():
            cls._dump_label(io, label)

    @classmethod
    def _load_coords(cls, io):
        cdim = io.load_int()
        assert -1 <= cdim <= 3
        pnts = io.load_cls(IS)
        dofs = io.load_cls(IS)
        crds = io.load_cls(Vec)
        if cdim > 0:
            crds.shape = (-1, cdim)
        return (pnts, dofs, crds)

    @classmethod
    def _dump_coords(cls, io, obj):
        pnts, dofs, crds = obj
        if len(crds.shape) == 2:
            cdim = crds.shape[-1]
        else:
            cdim = -1
        io.dump_int(cdim)
        io.dump_cls(IS, pnts)
        io.dump_cls(IS, dofs)
        io.dump_cls(Vec, crds)

    class Label(OrderedDict):
        name = None
        default = -1

    @classmethod
    def _load_label(cls, io):
        label = cls.Label()
        label.name = io.load_string(64)
        label.default = io.load_int()
        values = io.load_cls(IS)
        for value in values:
            points = io.load_cls(IS)
            label[value] = points
        return label

    @classmethod
    def _dump_label(cls, io, obj):
        label = obj
        io.dump_string(label.name, 64)
        io.dump_int(label.default)
        values = list(label.keys())
        values = np.asarray(values)
        io.dump_cls(IS, values)
        for points in label.values():
            io.dump_cls(IS, points)


@IO.register
class SSDC(object):

    CLASSID = 1211277

    def __init__(self):
        self.dm = None
        self.nsd = None
        self.coordpnts = None
        self.coorddofs = None
        self.coordinates = None
        self.bnd_label = None
        self.deg_label = None
        self.ref_label = None
        self.bnd_map = None
        self.deg_map = None
        self.ref_map = None
        self.deg = None
        self.dof = None

    @classmethod
    def load(cls, io):
        ssdc = cls()
        io.load_cid(cls)
        ssdc.dm = io.load_cls(DMPlex)
        ssdc.nsd = io.load_int()
        if ssdc.nsd > 0:
            ssdc.coordpnts = io.load_cls(IS)
            ssdc.coorddofs = io.load_cls(IS)
            ssdc.coordinates = io.load_cls(Vec)
            ssdc.coordinates.shape = (-1, ssdc.nsd)
        ssdc.bnd_label = cls._load_lbl(io)
        ssdc.deg_label = cls._load_lbl(io)
        ssdc.ref_label = cls._load_lbl(io)
        ssdc.bnd_map = cls._load_map(io)
        ssdc.deg_map = cls._load_map(io)
        ssdc.ref_map = cls._load_map(io)
        ssdc.deg = io.load_int()
        ssdc.dof = io.load_int()

        return ssdc

    @classmethod
    def dump(cls, io, obj):
        ssdc = obj
        io.dump_cid(cls)
        io.dump_cls(DMPlex, ssdc.dm)
        io.dump_int(ssdc.nsd)
        if ssdc.nsd > 0:
            io.dump_cls(IS, ssdc.coordpnts)
            io.dump_cls(IS, ssdc.coorddofs)
            io.dump_cls(Vec, ssdc.coordinates)
        cls._dump_lbl(io, ssdc.bnd_label)
        cls._dump_lbl(io, ssdc.deg_label)
        cls._dump_lbl(io, ssdc.ref_label)
        cls._dump_map(io, ssdc.bnd_map)
        cls._dump_map(io, ssdc.deg_map)
        cls._dump_map(io, ssdc.ref_map)
        io.dump_int(ssdc.deg)
        io.dump_int(ssdc.dof)

    @classmethod
    def _load_lbl(cls, io):
        return io.load_string(64) or None

    @classmethod
    def _dump_lbl(cls, io, obj):
        return io.dump_string(obj or '', 64)

    @classmethod
    def _load_map(cls, io):
        size = io.load_int()
        if size < 0:
            return None
        key = io.load_int(size)
        val = io.load_int(size)
        return OrderedDict(zip(key, val))

    @classmethod
    def _dump_map(cls, io, obj):
        if obj is None:
            io.dump_int(-1)
            return
        io.dump_int(len(obj))
        io.dump_int(obj.keys())
        io.dump_int(obj.values())


@IO.register
class Grid(object):

    CLASSID = 1211279

    @classmethod
    def load(cls, io):
        io.load_cid(cls)
        dim = io.load_size()
        degree = io.load_cls(IS)
        coords = io.load_cls(Vec)
        num = coords.size // dim
        coords.shape = (num, dim)
        assert np.sum((degree+1)**dim) == num
        return (degree, coords)

    @classmethod
    def dump(cls, io, obj):
        degree, coords = obj
        num, dim = np.shape(coords)
        assert np.sum((degree+1)**dim) == num
        io.dump_cid(cls)
        io.dump_size(dim)
        io.dump_cls(IS, degree)
        io.dump_cls(Vec, coords)


@IO.register
class Step(object):

    CLASSID = 1211278

    @classmethod
    def load(cls, io):
        io.load_cid(cls)
        step = io.load_int()
        time = io.load_real()
        data = io.load_cls(Vec)
        return (step, time, data)

    @classmethod
    def dump(cls, io, obj):
        step, time, data = obj
        io.dump_cid(cls)
        io.dump_int(step)
        io.dump_real(time)
        io.dump_cls(Vec, data)


@IO.register
class Avg(object):

    CLASSID = 1211280

    @classmethod
    def load(cls, io):
        io.load_cid(cls)
        init = io.load_real()
        time = io.load_real()
        data = io.load_cls(Vec)
        return (init, time, data)

    @classmethod
    def dump(cls, io, obj):
        init, time, data = obj
        io.dump_cid(cls)
        io.dump_real(init)
        io.dump_real(time)
        io.dump_cls(Vec, data)


def _main():
    # pylint: disable=import-outside-toplevel
    import sys

    indices = '32bit'
    precision = 'double'
    for arg in sys.argv[1:]:
        if arg in ('-i4', '-i32'):
            indices = '32bit'
            continue
        if arg in ('-i8', '-i64'):
            indices = '64bit'
            continue
        if arg in ('-f4', '-f32'):
            precision = 'single'
            continue
        if arg in ('-f8', '-f64'):
            precision = 'double'
            continue

        filename = arg
        io = IO(precision=precision, indices=indices)
        with io.open(filename, 'r'):
            cls = io.peek()
            while True:
                print("{0}: {1}".format(filename, cls.__name__))
                obj = io.load()
                del obj
                try:
                    cls = io.peek()
                except EOFError:
                    break


if __name__ == '__main__':
    _main()
