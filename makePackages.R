require(Rcpp)

cfiles = paste("Cpp/",list.files("Cpp/"),sep="")
Rfiles = paste("R/",list.files("R/"),sep="")
packageName="BayDR"
verNo=2.2
system(paste("rm -rf",packageName))
#  delete the original bayDR folder 
Rcpp.package.skeleton(name=packageName,cpp_files=cfiles,code_files=Rfiles,author="YL",
                      example_code=FALSE,email="yunlong.nie@stat.ubc.ca")


require(roxygen2)

roxygenise(package.dir="BayDR",overwrite=TRUE) # wirte some in the C++ files too, like in the DatDen_X_C

# add useDynLib(BayDR) at the end of NAMESPACE

system(paste("cp -r data ",paste0(packageName,"/")))
sink(file=paste0(packageName,"/NAMESPACE"),append=TRUE)
cat("useDynLib(",packageName,")",sep="")
sink()
file.copy(from="DESCRIPTION",to=paste0(packageName,"/"),overwrite=TRUE)
desp= readLines(con=paste0(packageName,"/","DESCRIPTION"))
desp[4] = paste("Version: ", verNo)
writeLines(desp,con=paste0(packageName,"/","DESCRIPTION"))
# then copy the description file and data folder, change the verion number in the description file. 

# command line to the the folder which bayDR locates

system(paste("R CMD build",packageName))
filenames = list.files("~/Dropbox/UBC/CinR/FileOnSever/Mypackages/Import/")
system(paste("R CMD check",filenames[grep(as.character(verNo),filenames)]))
system(paste("R CMD install",filenames[grep(as.character(verNo),filenames)]))



# R CMD build BayDR
# R CMD check BayDR
# R CMD INSTALL BayDR

