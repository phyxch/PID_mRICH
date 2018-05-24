#!/bin/csh
uname -a | cat > out.log
more /etc/redhat-release | cat >> out.log
lscpu | cat >> out.log
more /proc/meminfo | cat >> out.log

#use official installation on farm and ifarm
#setenv LD_LIBRARY_PATH /group/solid/solid_svn/evgen/cteq-pdf-1.0.4/lib
#setenv LIBRARY_PATH /group/solid/solid_svn/evgen/cteq-pdf-1.0.4/lib
source /group/eic/eic_svn/set_eic 1.3 
# source /group/eic/eic_svn/set_eic 2.2 

echo input list: $1 | cat >> out.log
echo input exe:  $2 | cat >> out.log

./$2 $1 >> out.log # database
# ./genMassHypo $1 >> out.log # database
# ./calLikelihood $1 >> out.log # likelihood


