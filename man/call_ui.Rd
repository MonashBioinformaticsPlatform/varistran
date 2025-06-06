% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_util.R
\name{call_ui}
\alias{call_ui}
\title{Call object with "request" if it is callable.}
\usage{
call_ui(ui, request)
}
\arguments{
\item{ui}{A UI object, or preferably a function(request) to produce a UI object.}

\item{request}{A request object to be passed to ui, if it is a function.}
}
\value{
A UI object.
}
\description{
This is used to support older Shiny UI code which doesn't wrap UI in \code{function(request) { ... }}.
}
