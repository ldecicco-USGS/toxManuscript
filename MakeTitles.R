
mktitl <- function(str){
  words = strsplit(str,' ')
  res = sapply(words[[1]],function(w){substring(w,1,1) = toupper(substring(w,1,1));w})
  paste(res,collapse=' ')
}

MakeTitles <- function(strvector) {
  strvector <- tolower(strvector)
  sapply(strvector,mktitl)
}