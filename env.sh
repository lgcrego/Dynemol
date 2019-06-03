#!/bin/bash
rm -r -f tmp_data 2> qdynamo.err
mkdir tmp_data
rm -r -f dos_trunk 2> qdynamo.err
mkdir dos_trunk
rm -r -f MO_trunk 2> qdynamo.err
mkdir MO_trunk
mkdir opt_trunk 2> qdynamo.err
mv opt_trunk/view_cost.dat opt_trunk/old_view_cost.dat 2> qdynamo.err
