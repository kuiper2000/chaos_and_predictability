#!bin/sh
git init
git add . 
git commit -m "Your message about the commit"
git remote add origin https://github.com/kuiper2000/chaos_and_predictability.io
git push -u origin gh-pages
git push origin gh-pages  
