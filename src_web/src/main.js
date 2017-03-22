import Vue from 'vue'
import App from './App.vue'
import Inputs from './Inputs.vue'
import Charts from './Charts.vue'
import Results from './Results.vue'

new Vue({
  el: '#app',
  render: h => h(App)
})

new Vue({
  el: '#inputFilesSection',
  render: h => h(Inputs)
})

new Vue({
  el: '#charts',
  render: h => h(Charts)
})

new Vue({
  el: '#resultsSelection',
  render: h => h(Results)
})

// TODO : to get a HighCharts component : https://github.com/vuejs/awesome-vue#libraries--plugins
