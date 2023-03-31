@echo off
echo Compiling...
comsolcompile SECM_standard.java
echo compilation complete.
echo simulating...
comsolbatch -inputfile SECM_standard.class -batchlog batlog.log
echo Process complete.
pause
