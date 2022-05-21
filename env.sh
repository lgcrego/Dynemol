#!/bin/bash

if [ $# -eq 0 ]
then 
	rm -r -f $DYNEMOLWORKDIR/dyn.trunk 2> $DYNEMOLDIR/qdynamo.err
	mkdir $DYNEMOLWORKDIR/dyn.trunk

	rm -r -f $DYNEMOLWORKDIR/dos.trunk 2> $DYNEMOLDIR/qdynamo.err
	mkdir $DYNEMOLWORKDIR/dos.trunk

	rm -r -f $DYNEMOLWORKDIR/MO.trunk 2> $DYNEMOLDIR/qdynamo.err
	mkdir $DYNEMOLWORKDIR/MO.trunk

	mkdir $DYNEMOLWORKDIR/opt.trunk 2> $DYNEMOLDIR/qdynamo.err
	mv $DYNEMOLWORKDIR/opt.trunk/view_cost.dat opt.trunk/old_view_cost.dat 2> $DYNEMOLDIR/qdynamo.err

	rm -r -f $DYNEMOLWORKDIR/log.trunk 2> $DYNEMOLDIR/qdynamo.err
	mkdir $DYNEMOLWORKDIR/log.trunk
fi

if [ "x$1" == 'xsave_cost_statement' ]
then
	paste $DYNEMOLWORKDIR/opt.trunk/view_cost.dat <(grep "eval(me)" $DYNEMOLWORKDIR/cost_tuning_EH.f | grep -v \!) > $DYNEMOLWORKDIR/opt.trunk/ga_cost.statement
fi
