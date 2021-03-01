const { createProxyMiddleware } = require("http-proxy-middleware");

// CRA has a proxyTimeout of 2 minutes, such that some long-running requests timeout with ERR_EMPTY_RESPONSE.
// https://create-react-app.dev/docs/proxying-api-requests-in-development/#configuring-the-proxy-manually
module.exports = function (app) {
  app.use(
    "/api",
    createProxyMiddleware({
      target: "http://api:5000",
      changeOrigin: true,
      // Set the timeout to 30min
      proxyTimeout: 1000 * 60 * 30,
      timeout: 1000 * 60 * 30,
    })
  );
};
