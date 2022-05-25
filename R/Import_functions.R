
#' @title Import Licor 6400 file
#' @description This functions allows to import the text file produced by LICOR as a data.frame
#' @param file File to import by the function
#' @param column_display The first lines of the file which are part of this list are displayed by this function after being imported.
#' @references Adapted from http://www.ericrscott.com/2018/01/17/li-cor-wrangling/
#' @return dataframe
#' @export
#'
#' @examples
f.import_licor6400<-function(file,column_display=c('Photo','Cond','PARi','Ci','Leaf_Barcode','Species','Tree Canopy','Age','file')){
  print(file)
  header_pattern <- "\"OPEN \\d\\.\\d\\.\\d"
  data_pattern <- "\\$STARTOFDATA\\$"
  text.raw <- readr::read_file(file)
  #splits into individual bouts
  raw_split <- stringr::str_split(text.raw, header_pattern, simplify = TRUE)

  #splits further to separate headers from actual data
  raw_split2 <- stringr::str_split(raw_split, data_pattern, simplify = FALSE)

  #extract just the second element, the actual data
  raw_split3 <- raw_split2 %>%
    map(`[`, 2) %>% #equivalent to doing raw_split2[[i]][2] for every element "i"
    flatten_chr() #converts to a vector

  #remove empty elements
  raw_split3 <- raw_split3[!is.na(raw_split3)]

  input <- raw_split3 %>%
    map(read_tsv, skip = 1)

  input.all <-do.call('rbind',lapply(X=input,FUN = function(x){x[as.vector(!is.na(x[,'HHMMSS'])),]}))#Supress comments
  imported=as.data.frame(input.all)
  imported=cbind(imported,file=rep(file,nrow(imported)))
  print(head(imported[,column_display]))
  return(imported)
}


#' @title Import Licor 6800 file
#' @description This functions allows to import the excel files produced by LICOR as a data.frame.
#' IMPORTANT: The excel files must be opened and saved before using this function (the Excel calculations are not done until the file is open, so the calculated colums will show 0s if not saved before being imported)
#' @param nskip_header Number of lines to skip in the Excel files to find the column names
#' @param nskip_data Number of lines to skip in the Excel files to find the data
#' @param file File path
#' @param column_display Column you want to display after the import to verufy if it worked correctly
#'
#' @examples
f.import_licor6800<-function(nskip_header=16,nskip_data=18,file,column_display=c('A','gsw','Qin','Ci','Species','Canopy','Pheno_Age','Barcode','file')){
  print(file)
  header=make.names(as.data.frame(readxl::read_excel(path = file,skip = nskip_header,n_max = 1,.name_repair = 'minimal',col_names = FALSE)))##'minimal' to speed up the import
  data_6800=as.data.frame(readxl::read_excel(path = file,skip = nskip_data,col_names = header,.name_repair = 'minimal'))
  data_6800[,'date']=data_6800[1,'date']
  data_6800=cbind(data_6800,file=rep(file,nrow(data_6800)))
  print(head(data_6800[,column_display]))
  return(data_6800)
}
