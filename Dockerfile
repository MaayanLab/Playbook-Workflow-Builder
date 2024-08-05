FROM node:21.3.0 as base
ENV PUPPETEER_SKIP_CHROMIUM_DOWNLOAD true
RUN echo "Installing system dependencies (git+puppeteer deps)..." \
  && apt-get update \
  && apt-get -y install \
    curl \
    git \
    gnupg \
  && curl --location --silent https://dl-ssl.google.com/linux/linux_signing_key.pub | apt-key add - \
  && sh -c 'echo "deb [arch=amd64] http://dl.google.com/linux/chrome/deb/ stable main" >> /etc/apt/sources.list.d/google.list' \
  && apt-get update \
  && apt-get install google-chrome-stable -y --no-install-recommends \
  && rm -rf /var/lib/apt/lists/*
RUN npm i -g ts-node
USER node
WORKDIR /app

FROM base as prepare_system
USER root
RUN echo "Installing system deps..." && apt-get -y update && apt-get -y install r-base python3-dev python3-pip python3-venv pkg-config libhdf5-dev && rm -rf /var/lib/apt/lists/*
ENV PYTHON_BIN="python3"
USER node

FROM prepare_system as prepare_r
USER root
COPY --chown=node:node cli/setup.R /app/setup.R
RUN echo "Running setup.R..." && R -e "source('/app/setup.R')" && rm /app/setup.R
USER node

FROM base as prepare_src
COPY --chown=node:node . /app

FROM prepare_src as prepare_package_json
RUN find /app -type d -name "node_modules" -exec rm -rf {} + \
  && find /app -type f -a \! \( -name "package.json" -o -name "package-lock.json" -o -name ".puppeteerrc.cjs" \) -delete \
  && find /app -type d -empty -delete

FROM prepare_src as prepare_requirements_txt
RUN find /app -type d -name "node_modules" -exec rm -rf {} + \
  && find /app -type f -a \! \( -name "requirements.txt" -o -name "enumerate-requirements.ts" \) -delete \
  && find /app -type d -empty -delete

FROM base as prepare_npm_i
COPY --from=prepare_package_json /app /app
RUN echo "Installing NodeJS dependencies..." && npm i

FROM prepare_npm_i as prepare_requirements_txt_complete
COPY --from=prepare_requirements_txt /app /app
RUN npm run codegen:requirements
USER root
RUN mv /app/requirements.txt /tmp/requirements.txt \
  && rm -r /app \
  && mkdir /app \
  && mv /tmp/requirements.txt /app \
  && chown node:node /app/requirements.txt
USER node

FROM prepare_src as prepare_build
COPY --from=prepare_npm_i /app /app
RUN echo "Building app..." && LANDING_PAGE=/graph/extend PUBLIC_URL=https://playbook-workflow-builder.cloud npm run build

FROM prepare_system as prepare_python
COPY --from=prepare_requirements_txt_complete /app /app
USER root
RUN echo "Installing python dependencies..." && python3 -m pip install --break-system-packages -r /app/requirements.txt && rm /app/requirements.txt
USER node

# TARGET: dev -- development environment with dependencies to run dev tools
FROM prepare_system as dev
USER root
RUN echo "Installing dev deps..." \
  && apt-get -y update \
  && apt-get -y install \
    librsvg2-bin \
    imagemagick \
    potrace \
    autoconf \
    pkg-config \
    build-essential \
    curl \
    libpng-dev \
  && ln -s /usr/bin/python3 /usr/bin/python \
  && curl -LO https://github.com/ImageMagick/ImageMagick/archive/refs/tags/7.1.1-26.tar.gz \
  && tar xzf 7.1.1-26.tar.gz \
  && rm 7.1.1-26.tar.gz \
  && sh ./ImageMagick-7.1.1-26/configure --prefix=/usr/local --with-bzlib=yes --with-fontconfig=yes --with-freetype=yes --with-gslib=yes --with-gvc=yes --with-jpeg=yes --with-jp2=yes --with-png=yes --with-tiff=yes --with-xml=yes --with-gs-font-dir=yes --with-rsvg=yes \
  && make -j && make install && ldconfig /usr/local/lib/ \
  && apt-get clean \
  && apt-get autoremove \
  && rm -rf /var/lib/apt/lists/*
USER node

COPY --from=amacneil/dbmate /usr/local/bin/dbmate /usr/local/bin/dbmate
ENV npm_config_cache=/app/.npm
COPY --from=prepare_r /usr/local/lib/ /usr/local/lib/
COPY --from=prepare_python /usr/local/lib/ /usr/local/lib/
CMD ["/bin/bash"]

# TARGET: app_minimal -- production server with dependencies to run just the webserver
FROM base as app_minimal
COPY --from=prepare_build /app /app
ENV PORT 3000
CMD ["npm", "run", "start"]

# TARGET: app -- production server with dependencies to run everything
FROM prepare_system as app
COPY --from=prepare_r /usr/local/lib/ /usr/local/lib/
COPY --from=prepare_python /usr/local/lib/ /usr/local/lib/
COPY --from=prepare_build /app /app
RUN chmod +x /app/cli/wes-worker.sh /app/cli/pwb.sh
ENV PORT 3000
ENV PORT 3005
CMD ["npm", "start"]
