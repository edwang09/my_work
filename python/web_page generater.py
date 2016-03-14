# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:08:41 2016

@author: mahone
"""
#import webbrowser
for i in range (2,101):
    f = open('page%d.html'%i,'w')
    
    message = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Untitled Document</title>
<link rel="stylesheet" type="text/css" href="pagestyle.css" />
</head>

<body>
<p>Introduction, pretest, and Profile</p>
<div class=boarder>
<div class=left1  title="apptosucceed"><img src="ats.png" alt="inavailable" width="330" height="70" /></div>
<div class=left2  title="Go to previous screen"><a  href=page%d.html><img src="left.png" alt="inavailable" width="70" height="70" /></a></div>
<div class=left3 title="Review topic"><img src="back.png" alt="inavailable" width="70" height="70" /></div>
<div class=left4 title="Learn more about this topic"><img src="learnmore.png" alt="inavailable" width="70" height="70" /></div>
<div class=left5 title="Take a break"><img src="chill.png" alt="inavailable" width="70" height="70" /></div>
<div class=left6 title="logout"><img src="logout.png" alt="inavailable" width="140" height="70" /></div>
<div class=left7 title="Need help? Ask The Wizard of Maze"><img src="ghost.png" alt="inavailable" width="80" height="75" /></div>
<div class=left8 title="Found the answer. Cancel my question"><img src="cancel.png" alt="inavailable" width="70" height="70" /></div>
<div class=left9 title="Skip this topic"><img src="skip.png" alt="inavailable" width="70" height="70" /></div>
<div class=left10 title="Submit answer or question"><img src="submit.png" alt="inavailable" width="70" height="70" /></div>
<div class=right title="Go to next screen"><a  href=page%d.html><img src="right.png" alt="inavailable" width="70" height="70"/></a></div>
</div>
<footer>
  <p>%d</p>
</footer>
</body>
</html>"""%(i-1,i+1,i)
    
    f.write(message)
    f.close()

#webbrowser.open_new_tab('helloworld.html')