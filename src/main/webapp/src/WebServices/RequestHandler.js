const RequestHandler = new function() {
  this.getRequest = function(url, /*, argumentToPass1, argumentToPass2, etc. */) {
    return new Promise(function (resolve, reject) {
      const req = new XMLHttpRequest();
      // Put the arguments after fcallback into a single arguments object
      req.arguments = Array.prototype.slice.call(arguments, 1);
      // Configure the behavior of the req object
      req.onload = function() {
        if (this.status >= 200 && this.status < 300) {
          resolve(this.responseText);
        } else {
          reject({
            status: this.status,
            statusText: this.statusText
          });
        }
      };
      req.onerror = function() {
        reject({
          status: this.status,
          statusText: xhr.statusText
        });
      };
      // req.open(method. url, async, user, password)
      req.open("GET", url, true);
      req.send(null);
    });
  };

  this.deleteRequest = function(url /*, argumentToPass1, argumentToPass2, etc. */) {
    return new Promise(function (resolve, reject) {
      const req = new XMLHttpRequest();
      // Put the arguments after fcallback into a single arguments object
      req.arguments = Array.prototype.slice.call(arguments, 1);
      // Configure the behavior of the req object
      req.onload = function() {
        if (this.status >= 200 && this.status < 300) {
          resolve(this.responseText);
        } else {
          reject({
            status: this.status,
            statusText: this.statusText
          });
        }
      };
      req.onerror = function() {
        reject({
          status: this.status,
          statusText: xhr.statusText
        });
      };
      // req.open(method. url, async, user, password)
      req.open("DELETE", url, true);
      req.send(null);
    });
  };

  this.postRequest = function(url, file,  /*, argumentToPass1, argumentToPass2, etc. */) {
    return new Promise(function (resolve, reject) {
      const req = new XMLHttpRequest();
      // Put the arguments after fcallback into a single arguments object
      req.arguments = Array.prototype.slice.call(arguments, 2);
      // Configure the behavior of the req object
      req.onload = function() {
        if (this.status >= 200 && this.status < 300) {
          resolve(this.responseText);
        } else {
          reject({
            status: this.status,
            statusText: this.statusText
          });
        }
      };
      req.onerror = function() {
        reject({
          status: this.status,
          statusText: xhr.statusText
        });
      };
      // req.open(method. url, async, user, password)
      req.open("POST", url, true);
      // Put the file into a FormData object needed to send the file through the request
      const fileToSend = new FormData();
      fileToSend.append("file", file);
      // Finally send the request containing the file
      req.send(fileToSend);
    });
  };
};

export default RequestHandler;
