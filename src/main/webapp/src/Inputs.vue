<template>
  <section id="inputFilesSection">
    <h1>Input files</h1>
    <label for="geolInput">Geologic.input : </label>
    <input type="file" name="geolInput" @change="onFileChange"id="geolInput" class="btn btn-default"/>
    <label for="morphoInput">Morphologic.input : </label>
    <input type="file" name="hydroInput" @change="onFileChange" id="morphoInput" class="btn btn-success"/>
    <label for="hydroInput">Hydrologic.input : </label>
    <input type="file" name="hydroInput" @change="onFileChange" id="hydroInput" class="btn btn-info"/>
    <label for="paramInput">Key_spatialized_param.param : </label>
    <input type="file" name="paramInput" @change="onFileChange" id="paramInput" class="btn btn-primary"/>
    <button class="btn btn-warning submit" @click="uploadFiles">Envoyer</button>
    <h3 id="inputFilesUploadSucces" v-if="allFilesUploaded">
      All files uploaded !
    </h3>
  </section>
</template>

<script>
import WS from './WebServices/WebServices.js'

export default {
  name: 'inputs',
  data () {
    return {
      geolInput: null,
      morphoInput: null,
      hydroInput: null,
      paramInput: null,
      allFilesUploaded: false,
    }
  },
  methods: {
    uploadFiles: function() {
      /* For each file input, verify if the user input a file. If so, upload it.
      After it : call this.fireMissingElementToSimulateRequest */
      if (this.geolInput != null) WS.postUploadGeolInput(this.geolInput).then(this.fireMissingElementToSimulateRequest).catch(function(err){console.error(err.statusText)});
      if (this.morphoInput != null) WS.postUploadMorphoInput(this.morphoInput).then(this.fireMissingElementToSimulateRequest).catch(function(err){console.error(err.statusText)});
      if (this.hydroInput != null) WS.postUploadHydroInput(this.hydroInput).then(this.fireMissingElementToSimulateRequest).catch(function(err){console.error(err.statusText)});
      if (this.paramInput != null) WS.postUploadParamInput(this.paramInput).then(this.fireMissingElementToSimulateRequest).catch(function(err){console.error(err.statusText)});
    },
    onFileChange: function(e){
      // Retrieve the file
      const file = e.target.files[0] || e.dataTransfer.files[0];
      // Retrieve the name of the input
      const modelData = e.srcElement.id;
      // Put the file into the data object
      this._data[modelData] = file;
    },
    // Callback for any uploadFileMethod
    fireMissingElementToSimulateRequest: function(){
      WS.getMissingElementToSimulate().then(this.updateAllFilesUploaded).catch(function(err){console.error(err.statusText)});
    },
    /* Called after receiving the responseCode from an upload file request.
    This method is used to update the allFilesUploaded boolean that display, or
    not, the "Successfully submitted !" message*/
    updateAllFilesUploaded: function(response) {
      const list = JSON.parse(response);
      (list.length == 0) ? this.allFilesUploaded = true : this.allFilesUploaded = false;
    }
  }
}
</script>

<style lang="scss">
  .submit {
    margin-top: 30px;
  }

  label {
    margin-top: 15px;
  }
</style>
