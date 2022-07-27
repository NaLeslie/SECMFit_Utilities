@echo off
echo Compiling...
comsolcompile SECM_Grid_Fit.java
echo Compilation complete.
echo Simulating...
comsolbatch -inputfile SECM_Grid_Fit.class -batchlog batlog.log
echo Process complete.
pause