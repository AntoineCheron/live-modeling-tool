import RequestHandler from './RequestHandler.js'

const API = new function () {
  // List of the webservices
  this.apiUrl = "API/";
  this.uploadGeolInput =  this.apiUrl + "uploadGeolInput";
  this.uploadHydroInput = this.apiUrl + "uploadHydroInput";
  this.uploadMorphoInput = this.apiUrl + "uploadMorphoInput";
  this.uploadParamInput = this.apiUrl + "uploadParamInput";
  this.simulate = this.apiUrl + "simulate";
  this.missingElementToSimulate = this.apiUrl + "missingElementToSimulate";
  this.result = this.apiUrl + "result/?name=";
  this.resultsFormat = this.apiUrl + "resultsFormat";

  this.postUploadGeolInput = function (file, callback) {
    RequestHandler.postRequest(this.uploadGeolInput, file, callback);
  };

  this.postUploadHydroInput = function (file, callback) {
    RequestHandler.postRequest(this.uploadHydroInput, file, callback);
  };

  this.postUploadMorphoInput = function (file, callback) {
    RequestHandler.postRequest(this.uploadMorphoInput, file, callback);
  };

  this.postUploadParamInput = function (file, callback) {
    RequestHandler.postRequest(this.uploadParamInput, file, callback);
  };

  this.getSimulate = function (callback) {
    RequestHandler.getRequest(this.simulate, callback);
  }

  this.getMissingElementToSimulate = function (callback) {
    RequestHandler.getRequest(this.missingElementToSimulate, callback);
  }

  this.getResult = function (name, callback) {
    RequestHandler.getRequest(`${this.simulate}${name}`, callback);
  }

  this.getResulsFormat = function (callback) {
    RequestHandler.getRequest(this.resultsFormat, callback);
  }

}

export default API;
