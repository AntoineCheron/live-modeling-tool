function Colors () {
  // Object containing all colors that will be used for the charts
  this.colors = [
    '#F44336',
    '#AB47BC',
    '#2196F3',
    '#26A69A',
    '#FFA726',
    '#A1887F',
    '#00E676',
    '#EC407A',
    '#7986CB',
    '#536DFE',
    '#26A69A',
    '#1DE9B6',
    '#FDD835',
    '#D4E157'
  ];
  this.nextIndex = -1;
}

Colors.prototype.nextColor = function () {
  // Make sure we never go out of the this.colors array
  this.nextIndex === this.colors.length-1 ? this.nextIndex=0 : this.nextIndex++;
  // Return the color
  return this.colors[this.nextIndex];
}

export default Colors;
