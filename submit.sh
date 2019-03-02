#!/bin/bash
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s13 --patch deep6 --array pa1 --overwrite --radial-fit-annulus 80 --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 15" --walltime 00:40:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s13 --patch deep5 --array pa1 --overwrite --radial-fit-annulus 80 --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 15" --walltime 00:40:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s13 --patch deep1 --array pa1 --overwrite --radial-fit-annulus 80 --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 15" --walltime 00:40:00 -t 80

mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch deep8 --array pa1 --overwrite --radial-fit-annulus 80 --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 15" --walltime 00:40:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch deep8 --array pa2 --overwrite --radial-fit-annulus 80 --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 15" --walltime 00:40:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch deep8 --array pa3 --overwrite --radial-fit-annulus 80 --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 15" --walltime 00:40:00 -t 80

mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s14 --patch deep56 --array pa1 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 5" --walltime 00:40:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s14 --patch deep56 --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 5" --walltime 00:40:00 -t 80

mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch deep56 --array pa1 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 5" --walltime 00:40:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch deep56 --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 5" --walltime 00:40:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch deep56 --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 200 --nsims 5" --walltime 00:40:00 -t 80


mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch boss --array pa1 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --nsims 5" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch boss --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --nsims 5" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s15 --patch boss --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --nsims 5" --walltime 01:00:00 -t 80



mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch000 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch001 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch002 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch003 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch004 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch005 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch006 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch007 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa2 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch008 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80


mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch000 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch001 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch002 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch003 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch004 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch005 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch006 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch007 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80
mpi_niagara 1 "python bin/make_covsqrt.py v1.3 act_mr3 --season s16 --patch cmb --array pa3 --overwrite --mask-version mr3c_20190215_pickupsub_190301 --mask-pad 400 --mask-patch patch008 --nsims 5 --dfact 4" --walltime 01:00:00 -t 80

