#!/bin/bash

#DYNEMOLWORKDIR="$(pwd)"
#export DYNEMOLWORKDIR

if [ $# -eq 0 ]
then 
	rm -r -f "$DYNEMOLWORKDIR"/dyn.trunk 2> "$DYNEMOLDIR"/qdynamo.err
	mkdir "$DYNEMOLWORKDIR"/dyn.trunk

	rm -r -f "$DYNEMOLWORKDIR"/dos.trunk 2> "$DYNEMOLDIR"/qdynamo.err
	mkdir "$DYNEMOLWORKDIR"/dos.trunk

	rm -r -f "$DYNEMOLWORKDIR"/MO.trunk 2> "$DYNEMOLDIR"/qdynamo.err
	mkdir "$DYNEMOLWORKDIR"/MO.trunk

	mkdir "$DYNEMOLWORKDIR"/opt.trunk 2> "$DYNEMOLDIR"/qdynamo.err
	mv "$DYNEMOLWORKDIR"/opt.trunk/ga_cost.statement opt.trunk/old_ga_cost.statement 2> "$DYNEMOLDIR"/qdynamo.err

	rm -r -f "$DYNEMOLWORKDIR"/log.trunk 2> "$DYNEMOLDIR"/qdynamo.err
	mkdir "$DYNEMOLWORKDIR"/log.trunk

	rm -r -f "$DYNEMOLWORKDIR"/ancillary.trunk 2> "$DYNEMOLDIR"/qdynamo.err
	mkdir "$DYNEMOLWORKDIR"/ancillary.trunk
	mkdir "$DYNEMOLWORKDIR"/ancillary.trunk/configs
fi

if [ "x$1" == 'xsave_cost_statement' ]
then
	paste "$DYNEMOLWORKDIR"/opt.trunk/view_cost.dat <(grep "eval(me)" "$DYNEMOLDIR"/cost_tuning_EH.f | grep -v \!) > "$DYNEMOLWORKDIR"/opt.trunk/ga_cost.statement
    rm "$DYNEMOLWORKDIR"/opt.trunk/view_cost.dat
fi
