// Popup scripts for survey figures and tables

function popup(mylink, windowname) 
{                                 
   if (! window.focus)return true; 
      var href; 
      if (typeof(mylink) == 'string') 
         href=mylink; 
      else
         href=mylink.href; 
         window.open(href, windowname, 'toolbar=1,statusbar=1,menubar=1,scrollbars=1,location=1,width=900,height=900,left=20,right=20,resizable=1'); 
         return false; 
}

function popup640(mylink, windowname) 
{                                 
   if (! window.focus)return true; 
      var href; 
      if (typeof(mylink) == 'string') 
         href=mylink; 
      else
         href=mylink.href; 
         window.open(href, windowname, 'toolbar=1,statusbar=1,menubar=1,scrollbars=1,location=1,width=640,height=640,left=20,right=20,resizable=1'); 
         return false; 
}

function popup640x480(mylink, windowname) 
{                                 
   if (! window.focus)return true; 
      var href; 
      if (typeof(mylink) == 'string') 
         href=mylink; 
      else
         href=mylink.href; 
         window.open(href, windowname, 'toolbar=1,statusbar=1,menubar=1,scrollbars=1,location=1,width=640,height=480,left=20,right=20,resizable=1'); 
         return false; 
}

function popup400(mylink, windowname) 
{                                 
   if (! window.focus)return true; 
      var href; 
      if (typeof(mylink) == 'string') 
         href=mylink; 
      else
         href=mylink.href; 
         window.open(href, windowname, 'toolbar=1,statusbar=1,menubar=1,scrollbars=1,location=1,width=380,height=800,left=20,right=20,resizable=1'); 
         return false; 
}

function popupWide(mylink, windowname)
{
   if (! window.focus)return true;
      var href;
      if (typeof(mylink) == 'string')
         href=mylink;
      else
         href=mylink.href;
         window.open(href, windowname, 'toolbar=1,statusbar=1,menubar=1,scrollbars=1,location=1,width=1020,height=900,left=20,right=20,resizable=1');
         return false;
}

