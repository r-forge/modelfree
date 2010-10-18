// JavaScript Document

// Be default tracker number is for David's website
var tracknumber = "UA-10105990-1";
var kwebsite = /zychaluk/;

// Is this instead Kamila's website? Then, use Kamila's tracker number
var matchPos = location.href.search(kwebsite);
if(matchPos != -1) tracknumber = "UA-10105990-2"; 


 try {
 	var pageTracker = _gat._getTracker(tracknumber);
 	pageTracker._trackPageview();
 } catch(err) {}