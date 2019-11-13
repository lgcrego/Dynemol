#!/bin/bash

if [ $# -eq 0 ]
then 
	rm -r -f dyn_trunk 2> qdynamo.err
	mkdir dyn_trunk

	rm -r -f dos_trunk 2> qdynamo.err
	mkdir dos_trunk

	rm -r -f MO_trunk 2> qdynamo.err

	mkdir MO_trunk

	mkdir opt_trunk 2> qdynamo.err
	mv opt_trunk/view_cost.dat opt_trunk/old_view_cost.dat 2> qdynamo.err
fi

var = $1
if [ "x$1" == 'xsave_cost_statement' ]
then
	paste opt_trunk/view_cost.dat <(grep "eval(me)" cost_tuning_EH.f | grep -v \!) > opt_trunk/ga_cost.statement
fi
