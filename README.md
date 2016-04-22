jsMetal
=======

Distance alignment metrics in JavaScript

Deployment
==========

The javascript files in scripts/ need to be combined into a single file. You could minify them,
or else you could just do:

        rm -f scripts/webmetal.js && cat scripts/*.js > scripts/webmetal.js

Then serve up this directory with any webserver and you should be set.


