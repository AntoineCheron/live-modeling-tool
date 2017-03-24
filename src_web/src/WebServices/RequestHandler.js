const RequestHandler = new function() {
  this.getRequest = function(url, fcallback /*, argumentToPass1, argumentToPass2, etc. */) {
    const req = new XMLHttpRequest();
    // Put the arguments after fcallback into a single arguments object
    req.arguments = Array.prototype.slice.call(arguments, 2);
    // Configure the behavior of the req object
    req.callback = fcallback;
    req.onload = this.xhrSuccess;
    req.onerror = this.xhrError;
    // req.open(method. url, async, user, password)
    req.open("GET", url, true);
    req.send(null);
  };

  this.deleteRequest = function(url, fcallback /*, argumentToPass1, argumentToPass2, etc. */) {
    const req = new XMLHttpRequest();
    // Put the arguments after fcallback into a single arguments object
    req.arguments = Array.prototype.slice.call(arguments, 2);
    // Configure the behavior of the req object
    req.callback = fcallback;
    req.onload = this.xhrSuccess;
    req.onerror = this.xhrError;
    // req.open(method. url, async, user, password)
    req.open("DELETE", url, true);
    req.send(null);
  };

  this.postRequest = function(url, file, fcallback,  /*, argumentToPass1, argumentToPass2, etc. */) {
    const req = new XMLHttpRequest();
    // Put the arguments after fcallback into a single arguments object
    req.arguments = Array.prototype.slice.call(arguments, 3);
    // Configure the behavior of the req object
    req.callback = fcallback;
    req.onload = this.xhrSuccess;
    req.onerror = this.xhrError;
    // req.open(method. url, async, user, password)
    req.open("POST", url, true);
    // Put the file into a FormData object needed to send the file through the request
    const fileToSend = new FormData();
    fileToSend.append("file", file);
    // Finally send the request containing the file
    req.send(fileToSend);
  }

  this.xhrSuccess = function() {
    // In this case this is the req object
    this.callback(this);
  };

  this.xhrError = function() {
    // In this case this is the req object
    console.error(this.statusText);
  };
};

export default RequestHandler;
