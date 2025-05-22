library(dplyr)
library(stringr)
if (FALSE){
  dataFile = "data/test_data/test_compound_classes_4.csv"
  samplesFile = "data/test_data/test_samples_op1.csv"
  first_row_groups = NULL
  path_to_save = NULL
  PeakArea_columns_pattern = "Norm..Area.."
  samples_order_name = F
  Tags = F
}

##############
prepare_metabo_tables <- function(dataFile,
                                  samplesFile,
                                  first_row_groups = NULL,
                                  path_for_saving_files = NULL,
                                  PeakArea_columns_pattern = "Norm..Area..",
                                  samples_order_name = F,
                                  Tags = F) {
  
  
  # load LCMS data
  # fix the first column name
  if (is.data.frame(dataFile)) {
    lc.raw.data <- dataFile
  } else {
    lc.raw.data <- read.csv(dataFile, encoding = "UTF8")
    #names(lc.raw.data) = gsub("?..","",names(lc.raw.data))
  }
  lc.data <- lc.raw.data
  if (Tags){
    lc.data <- lc.data %>%
      filter(Tags!="")
  }

  lc.data$Name <- iconv(lc.data$Name, "latin1", "UTF-8")
  
  ### if the column samples names *did change manually*
  if (is.null(nrow(samplesFile))&length(samplesFile)!=1) {
    
    
    lc.for.met<-lc.data%>%
      select(all_of(c("Comp_Class","Name", samplesFile)))%>%
      mutate(Name=sub("std", "STD", Name))
    
  } else {
    ### if the column samples names *didnt change manually*
    # load samples list
    if (is.data.frame(samplesFile)) {
      samples <- samplesFile
    } else {
      samples <- read.csv(samplesFile, encoding = "UTF8")
      names(samples)[1] = "x"
      #names(samples) = gsub("?..","",names(samples))
    }
    
    #############samples$group = c(rep("M1",3), rep("M3",3), rep("M4",3), rep("M5",3), rep("EV",4))
    
    # focus only on "Norm.Area" fields of samples:
    # >> extract the sample id from column name
    samples_pattern = paste(PeakArea_columns_pattern,samples$x, sep = "", collapse = "|")
    
    LCMS.NormArea.orig<-grep(samples_pattern ,names(lc.data), value=TRUE, ignore.case = TRUE)
    samplesOrder<-samples$x[match(sub("Norm\\.\\.Area\\.\\.(.*)\\.raw\\.\\..*", "\\1", LCMS.NormArea.orig),
                                       samples[,1])]
    
    # >> arrange samples rows by the order of samples in LCMS
    samples.ordered.df<-samples[match(samplesOrder,samples$x),]
    samplesNames<-samples.ordered.df$sample
    
    #select only column of metabolite's name and Norm Area of samples
    # and rearrange column order by samplesOrder
    # and replace LCMS Norm area names with samples names from samples file
    # and replcae std text with upper case [sivan: ask Yonathan about it]
    lc.for.met<-lc.data%>%
      select(all_of(c("Comp_Class","Name", LCMS.NormArea.orig)))%>%
      rename_at(vars(LCMS.NormArea.orig), ~samplesNames)%>%
      mutate(Name=sub("std", "STD", Name))
    
  }
  
  
  # if you want to get only <samples_order_name>
  if (samples_order_name == T) {
    print(names(lc.for.met)[-1])
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  
  
  # make non duplicates names and classes
  non_duplicates <- function(string) {
    mstring <- make.unique(as.character(string), sep=" X" )
    tmp <- !duplicated(string)
    for (i in 1:length(mstring[tmp])){
      mstring[tmp][i]<-ifelse(string[tmp][i] %in% string[duplicated(string)]
                              , gsub("(.*)","\\1 X0", mstring[tmp][i])
                              , mstring[tmp][i]
      )
    }
    end <- sub(".* X([0-9]+)","\\1",grep(" X([0-9]*)$",mstring,value=T) )
    beg <- sub("(.* X)[0-9]+","\\1",grep(" X([0-9]*)$",mstring,value=T) )
    newend <- as.numeric(end)+1
    mstring[grep(" X([0-9]*)$",mstring)] <- paste0(beg,newend)
    mstring 
  }
  lc.for.met$Name = non_duplicates(lc.for.met$Name)
  lc.for.met$Comp_Class = non_duplicates(lc.for.met$Comp_Class)
  
  
  # add first_row_groups to the first row
  if (is.null(first_row_groups)) {
    first_row_groups = samples.ordered.df$group
  }
  lc.for.met = rbind(
    c("","",first_row_groups),
    lc.for.met
  )
  
  # if there is no path to save, the function will return all tables as vectors
  if (is.null(path_for_saving_files) == F) {
    # create new directory in path_for_saving_files and save all tables
    # remove  "/" if needed
    if (str_sub(path_for_saving_files, nchar(path_for_saving_files)) == "/") {
      path_for_saving_files = str_sub(path_for_saving_files, 1,nchar(path_for_saving_files)-1)
    }
    
    # create new directory for results
    date.suffix = gsub(" ","",format(Sys.time(), "%d %m %y"))
    path.for.output.folder = paste0(path_for_saving_files,"/files_for_metaboanalyst_",date.suffix)
    dir.create(path.for.output.folder)
    
    # save class file
    write.csv(lc.for.met[-2],
              sprintf("%s/comp.class.for.metabo.%s.csv",path.for.output.folder, date.suffix),
              row.names = F)
    # save name file
    write.csv(lc.for.met[,-1],
              sprintf("%s/comp.name.for.metabo.%s.csv",path.for.output.folder, date.suffix),
              row.names = F)
    # save name2class file
    write.csv(lc.for.met[-1,1:2],
              sprintf("%s/name2class.table.%s.csv",path.for.output.folder, date.suffix),
              row.names = F)
    
    message(sprintf("\n\n   *  files saved in '%s' folder  *\n\n",path.for.output.folder))
  }
  return(list(comp_name = lc.for.met[,-1],
              comp_class = lc.for.met[,-2],
              name2class = lc.for.met[-1,1:2]))
}

