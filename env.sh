#!/bin/bash

if [ $# -eq 0 ]
then 
	rm -r -f dyn.trunk 2> qdynamo.err
	mkdir dyn.trunk

	rm -r -f dos.trunk 2> qdynamo.err
	mkdir dos.trunk

	rm -r -f MO.trunk 2> qdynamo.err
	mkdir MO.trunk

	mkdir opt.trunk 2> qdynamo.err
	mv opt.trunk/view_cost.dat opt.trunk/old_view_cost.dat 2> qdynamo.err
fi

if [ "x$1" == 'xsave_cost_statement' ]
then
	paste opt.trunk/view_cost.dat <(grep "eval(me)" cost_tuning_EH.f | grep -v \!) > opt.trunk/ga_cost.statement
fi
