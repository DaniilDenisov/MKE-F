FROM mcr.microsoft.com/playwright:v1.64.0-noble
RUN npm install --prefix /opt/mkef-tests playwright@1.64.0
ENV MKEF_PLAYWRIGHT_MODULE=/opt/mkef-tests/node_modules/playwright
WORKDIR /app
CMD ["node", "scripts/test_sections_browser.cjs"]
