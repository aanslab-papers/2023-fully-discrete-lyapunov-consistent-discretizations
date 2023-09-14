PARAMS = -Rd 2 -nu 0.15
TS3BS = -rrk 1
BISECT = -rrk_solver 0
TSADAPTMAX = -ts_adapt_dt_max 1e-10
EXEC = mpiexec

srun:
	@export EXEC=srun
.PHONY: srun


run: data1 data2 data3 data4 data5 data6 data7 data8 data9 data10 data11 data12 data13

data/pde-scaled-si-eq.txt: model-pde-scaled-si.bin
	-@mkdir -p data
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin -Rd 1.5 -nu 1.1 -example 3 -ts_max_steps 0 -eq > $@
data1: data/pde-scaled-si-eq.txt

data/pde-scaled-si.txt: model-pde-scaled-si.bin
	-@mkdir -p data
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin -Rd 1.5 -nu 1.1 -example 3 -ts_max_time 27 $(TS3BS) $(BISECT) -max_cases_monitor > $@
data2: data/pde-scaled-si.txt

data/pde-scaled-si-rrk.txt: model-pde-scaled-si.bin
	-@mkdir -p data
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin -Rd 1.5 -nu 1.1 -example 3 -ts_max_time 27 $(TS3BS) $(BISECT) -rrk_monitor > $@
data3: data/pde-scaled-si-rrk.txt

data/figure_si/1/SCALED-SI.pvd: model-pde-scaled-si.bin
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin $(PARAMS) $(TS3BS) $(BISECT) -example 1 -vtk 1
	-@mkdir -p data/figure_si/1/
	-@mv SCALED-SI.pvd* data/figure_si/1/
data4: data/figure_si/1/SCALED-SI.pvd

data/figure_si/2/SCALED-SI.pvd: model-pde-scaled-si.bin
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin $(PARAMS) $(TS3BS) $(BISECT) -example 2 -vtk 1
	-@mkdir -p data/figure_si/2/
	-@mv SCALED-SI.pvd* data/figure_si/2/
data5: data/figure_si/2/SCALED-SI.pvd

data/figure_si/3/SCALED-SI.pvd: model-pde-scaled-si.bin
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin $(PARAMS) $(TS3BS) $(BISECT) -example 3 -vtk 1
	-@mkdir -p data/figure_si/3/
	-@mv SCALED-SI.pvd* data/figure_si/3/
data6: data/figure_si/3/SCALED-SI.pvd

data/figure_si/4/SCALED-SI.pvd: model-pde-scaled-si.bin
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin $(PARAMS) $(TS3BS) $(BISECT) -example 4 -vtk 1
	-@mkdir -p data/figure_si/4/
	-@mv SCALED-SI.pvd* data/figure_si/4/
data7: data/figure_si/4/SCALED-SI.pvd

data/figure_si/5/SCALED-SI.pvd: model-pde-scaled-si.bin
	-@$(call $(EXEC),1) ./model-pde-scaled-si.bin $(PARAMS) $(TS3BS) $(BISECT) -example 5 -vtk 1
	-@mkdir -p data/figure_si/5/
	-@mv SCALED-SI.pvd* data/figure_si/5/
data8: data/figure_si/5/SCALED-SI.pvd

data/pde-sei-eq.txt: model-pde-sei.bin
	-@mkdir -p data
	-@$(call $(EXEC),1) ./model-pde-sei.bin -example 1 -ts_max_steps 0 -eq > $@
data9: data/pde-sei-eq.txt

data/pde-sei.txt: model-pde-sei.bin
	-@mkdir -p data
	-@$(call $(EXEC),1) ./model-pde-sei.bin -example 1 $(TS3BS) $(BISECT) -max_cases_monitor > $@
data10: data/pde-sei.txt

data/pde-sei-rrk.txt: model-pde-sei.bin
	-@mkdir -p data
	-@$(call $(EXEC),1) ./model-pde-sei.bin -example 1 $(TS3BS) $(BISECT) -rrk_monitor > $@
data11: data/pde-sei-rrk.txt

data/figure_sei/sei/SEI.pvd: model-pde-sei.bin
	-@$(call $(EXEC),1) ./model-pde-sei.bin -example 1 $(TS3BS) $(BISECT) -ts_max_time 2 -vtk 1
	-@mkdir -p data/figure_sei/sei/
	-@mv SEI.pvd* data/figure_sei/sei/
data12: data/figure_sei/sei/SEI.pvd

data/m1virus-1D-rrk.txt: model-pde-m1virus.bin
	-@mkdir -p data
	-@$(call $(EXEC),1) ./model-pde-m1virus.bin -nel 32 -deg 2 -case 1 -ts_max_time 400 $(TS3BS) $(BISECT) $(TSADAPTMAX) -rrk_monitor > $@
data13: data/m1virus-1D-rrk.txt


