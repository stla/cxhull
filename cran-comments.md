## CRAN checks

The CRAN checks results find an error with the old release of R. That's 
because of an example which calls `grDevices::hcl.colors` with the palette 
named `Rocket`, and this palette was not available in R-4.0. So I modified 
this example by including this code:

```
if(getRversion() < "4.1.0"){
  palette <- "Viridis"
}else{
  palette <- "Rocket"
}
```


## Test environments

* Windows 10, R-4.1.2 
* win-builder (devel)
* mac-builder
* Ubuntu 20, via Github action


## R CMD check results

Status: OK


## Reverse dependencies

There's no consequence on the reverse dependencies: this new version only 
provides a new function.

