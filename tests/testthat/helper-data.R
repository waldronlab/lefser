load_lefser_dataset <- function(name) {
    dataenv <- new.env(parent = emptyenv())
    utils::data(list = name, package = "lefser", envir = dataenv)
    dataenv[[name]]
}
