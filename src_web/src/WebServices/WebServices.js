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

  this.postUploadGeolInput = function (file) {
    return RequestHandler.postRequest(this.uploadGeolInput, file);
  };

  this.postUploadHydroInput = function (file) {
    return RequestHandler.postRequest(this.uploadHydroInput, file);
  };

  this.postUploadMorphoInput = function (file) {
    return RequestHandler.postRequest(this.uploadMorphoInput, file);
  };

  this.postUploadParamInput = function (file) {
    return RequestHandler.postRequest(this.uploadParamInput, file);
  };

  this.getSimulate = function () {
    return RequestHandler.getRequest(this.simulate);
  }

  this.getMissingElementToSimulate = function () {
    return RequestHandler.getRequest(this.missingElementToSimulate);
  }

  this.getResult = function (name) {
    return RequestHandler.getRequest(`${this.result}${name}`);
  }

  this.getResulsFormat = function () {
    return RequestHandler.getRequest(this.resultsFormat);
  }

}

export default API;
