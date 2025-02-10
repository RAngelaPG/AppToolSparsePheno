.libPaths("./R-Portable/App/R-Portable/library")
message('library paths:\n', paste('... ', .libPaths(), sep='', collapse='\n'))
#chrome.portable = file.path(getwd(),'GoogleChromePortable/App/Chrome-bin/chrome.exe')
#launch.browser = function(appUrl, browser.path=chrome.portable) {
#    message('Browser path: ', browser.path)
#    shell(sprintf('"%s" --app=%s', browser.path, appUrl))
#}
shiny::runApp("./shiny/",port=8886,launch.browser=TRUE)
#shiny::runApp("./shiny/",launch.browser=launch.browser)
