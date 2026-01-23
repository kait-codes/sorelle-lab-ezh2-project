# Initial package install:
## Source - https://stackoverflow.com/a
## Posted by Matthew
## Retrieved 2026-01-12, License - CC BY-SA 3.0

using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

bioc_using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    BiocManager::install(need)
    lapply(need,require,character.only=TRUE)
  }
}