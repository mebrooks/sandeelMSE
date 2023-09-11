stocks = c("SAN_area_1r_C", 
					 "SAN_area_2r_C", 
					 "SAN_area_3r_C", 
					 "SAN_area_4_C")

for(stock.name in stocks){
	wd = file.path("assessments", stock.name)
	source(file.path(wd, "run_assessment.R"))
}
