
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

#' @title Import Licor 6800 excel file
#' @description This functions allows to import the excel file produced by LICOR as a data.frame. The files have to be open in Excel and saved before using his function so the result of the formula are calculated. The formula are sotred into the Excel file but not computed until the file is open.
#'
#' @inherit f.import_licor6400
#'
#' @return dataframe
#' @export
#'
#' @examples
f.import_licor6800<-function(file,column_display=c('A','gsw','Qin','Ci','Species','Canopy','Pheno_Age','Barcode','file')){
  print(file)
  header=make.names(as.data.frame(readxl::read_excel(path = file,skip = 16,n_max = 1,col_names = FALSE)))
  data_6800=as.data.frame(readxl::read_excel(path = file,skip = 18,col_names = header))
  data_6800[,'date']=data_6800[1,'date']
  data_6800=cbind(data_6800,file=rep(file,nrow(data_6800)))
  print(head(data_6800[,column_display]))
  return(data_6800)
}
