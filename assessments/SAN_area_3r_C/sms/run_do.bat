M:
cd "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\"
del /Q mdDir\*.*
del /Q retro\*.*
del /Q mcmc*.*
del /Q *.?0?
del /Q *.out
del /Q *.mcm
del /Q *.bin
del /Q admodel.*
del /Q *.csv
del /Q *.std
del /Q *.bar
del /Q *.mc2
del /Q *.cor
del /Q *.psv
del /Q *.ecm
del /Q *.xls
del /Q *.html
del /Q mcout*.all
del /Q *.wmf
del /Q *.png
del /Q *.lg
del /Q *.log
del /Q ud.dat
del /Q HCR_prob.dat
del /Q HCR_yield.dat
del /Q HCR_SSB.dat
del /Q *.par
del /Q *.rep
del /Q *.hst
del /Q *.eva
del /Q *.tmp
M:
cd "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\"
del /f "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\sms.rep 
del /f "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\sms.par 
del /f "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\sms.std 
sms  -nox -ind run_ms0.dat 
copy /Y "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\sms.par" "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\run_ms0.par" 
copy /Y "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\sms.rep" "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\run_ms0.rep" 
copy /Y "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\sms.std" "M:\Tobis\Tobis_assessment\SMS_2022\SAN-area-3r\run_ms0.std" 
