import plotComponentFactory from 'react-plotly.js/factory';
// Only use Plotly Geo to reduce bundle size significantly
// @ts-ignore
import Plotly from 'plotly.js-dist';
// Solution for using plotly.js-dist with react-plotly.js: https://github.com/plotly/react-plotly.js/issues/143
export const PlotComponent = plotComponentFactory(Plotly);
