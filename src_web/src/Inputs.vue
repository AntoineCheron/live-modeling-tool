<template>
  <section id="inputFilesSection">
    <div class="col-lg-3">
      <h1>Input files</h1>
      <label for="geolInput">Geologic.input : </label>
      <input type="file" name="geolInput" @change="onFileChange"id="geolInput" class="btn btn-default"/>
      <label for="morphoInput">Morphologic.input : </label>
      <input type="file" name="hydroInput" v-model="morphoInput" id="morphoInput" class="btn btn-success"/>
      <label for="hydroInput">Hydrologic.input : </label>
      <input type="file" name="hydroInput" v-model="hydroInput" id="hydroInput" class="btn btn-info"/>
      <label for="paramInput">Key_spatialized_param.param : </label>
      <input type="file" name="paramInput" v-model="paramInput" id="paramInput" class="btn btn-primary"/>
      <button class="btn btn-warning submit" v-on:click="uploadFiles">Envoyer</button>
      <p id="inputFilesUploadSucces" v-if="allFilesUploaded">
        Successfully submitted !
      </p>
    </div>
  </section>
</template>

<script>
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
      alert('hello');
    },
    onFileChange: function(e){
      // Retrieve the file
      const file = e.target.files[0] || e.dataTransfer.files[0];
      // Retrieve the name of the input
      const modelData = e.srcElement.id;
      // Put the file into the data object
      this._data[modelData] = file;
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
