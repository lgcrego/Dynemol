#!/bin/bash
rm -r -f tmp_data 2> qdynamo.err
mkdir tmp_data
rm -r -f DOS_trunk 2> qdynamo.err
mkdir dos_trunk
mkdir opt_trunk 2> qdynamo.err
mv OPT_trunk/view_cost.dat OPT_trunk/old_view_cost.dat 2> qdynamo.err
