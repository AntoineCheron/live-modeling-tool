<template>
  <section id="inputFilesSection">
    <div class="col-lg-3">
      <h1>Input files</h1>
      <form name="inputFiles" id="inputFiles" method="post" enctype="multipart/form-data" action="#">
        <label for="geolInput">Geologic.input : </label>
        <input type="file" name="geolInput" id="geolInput" class="btn btn-default"/>
        <label for="morphoInput">Morphologic.input : </label>
        <input type="file" name="hydroInput" id="morphoInput" class="btn btn-success"/>
        <label for="hydroInput">Hydrologic.input : </label>
        <input type="file" name="hydroInput" id="hydroInput" class="btn btn-info"/>
        <label for="paramInput">Key_spatialized_param.param : </label>
        <input type="file" name="paramInput" id="paramInput" class="btn btn-primary"/>
        <input type="submit" class="btn btn-warning"/>
      </form>
      <p id="inputFilesUploadSucces" class="hide">
        Successfully submitted
      </p>
    </div>
  </section>
</template>

<script>
export default {
  name: 'inputs',
  data () {
    return {
      geologic: 'geol',
      morphologic: 'morpho',
      hydrologic: 'hydro',
      param: 'key_spatialized_param',
    }
  }
}
</script>

<style lang="scss">
  #inputFiles .btn-warning {
    margin-top: 30px;
  }
</style>
