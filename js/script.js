/* Author: Matt Settles

*/
// Scroll to location
$(document).ready(function(){
  // This code is executed after the DOM has been completely loaded
  $('nav a,footer a.up').click(function(e){
  // If a link has been clicked, scroll the page to the link's hash target:
    $.scrollTo( this.hash || 0, 1500);
    e.preventDefault();
  });
});

// Google Maps
function gmaps.initialize() {
  var latlng = new google.maps.LatLng(-34.397, 150.644);
  var myOptions = {
    zoom: 8,
    center: latlng,
    mapTypeId: google.maps.MapTypeId.ROADMAP
  };
  var map = new google.maps.Map(document.getElementById("map_canvas"),
      myOptions);
}










