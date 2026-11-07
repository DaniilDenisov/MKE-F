FROM mcr.microsoft.com/playwright:v1.64.0-noble AS browser
RUN npm install --prefix /opt/mkef-tests playwright@1.64.0
ENV MKEF_PLAYWRIGHT_MODULE=/opt/mkef-tests/node_modules/playwright
WORKDIR /app
CMD ["node", "tests/browser/test_sections_browser.cjs"]

# Release acceptance also verifies the committed bundle against package-lock.json.
FROM browser AS release
COPY package.json package-lock.json /opt/mkef-assets/
RUN npm ci --prefix /opt/mkef-assets
ENV NODE_PATH=/opt/mkef-assets/node_modules
