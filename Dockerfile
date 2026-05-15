FROM maayanlab/base:1.4.0 AS base
USER root
ENV PUPPETEER_SKIP_CHROMIUM_DOWNLOAD=true
RUN echo "Installing system dependencies (git+puppeteer deps)..." \
  && apt-get update \
  && apt-get -y install \
    curl \
    git \
    gnupg \
    texlive-latex-extra \
    texlive-luatex \
    texlive-bibtex-extra \
    latexmk \
    biber \
  && tlmgr init-usertree \
  && tlmgr --usermode --repository https://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2023/tlnet-final install sourcesanspro \
  && curl --location --silent https://dl-ssl.google.com/linux/linux_signing_key.pub | apt-key add - \
  && sh -c 'echo "deb [arch=amd64] http://dl.google.com/linux/chrome/deb/ stable main" >> /etc/apt/sources.list.d/google.list' \
  && apt-get update \
  && apt-get install google-chrome-stable -y --no-install-recommends \
  && rm -rf /var/lib/apt/lists/*
ENV PUPPETEER_EXECUTABLE_PATH="/usr/bin/google-chrome"
USER ubuntu
RUN NODE_VERSION=21.3 /install.sh
RUN npm i -g ts-node

FROM base AS prepare_src
COPY --chown=ubuntu . app

FROM prepare_src AS prepare_package_json
RUN find app -type d -name "node_modules" -exec rm -rf {} + \
  && find app -type f -a \! \( -name "package.json" -o -name "package-lock.json" -o -name ".puppeteerrc.cjs" \) -delete \
  && find app -type d -empty -delete

FROM prepare_src AS prepare_requirements_txt
RUN find app -type d -name "node_modules" -exec rm -rf {} + \
  && find app -type f -a \! \( -name "requirements.txt" -o -name "enumerate-requirements.ts" \) -delete \
  && find app -type d -empty -delete

FROM base AS prepare_npm_i
COPY --from=prepare_package_json --chown=ubuntu /home/ubuntu/app app
RUN echo "Installing NodeJS dependencies..." && cd app && npm i

FROM prepare_npm_i AS prepare_requirements_txt_complete
COPY --from=prepare_requirements_txt --chown=ubuntu /home/ubuntu/app app
RUN cd app && npm run codegen:requirements
RUN mv app/requirements.txt /tmp/requirements.txt \
  && rm -r app \
  && mkdir app \
  && mv /tmp/requirements.txt app \
  && chown ubuntu app/requirements.txt

FROM prepare_src AS prepare_build
COPY --from=prepare_npm_i --chown=ubuntu /home/ubuntu/app app
RUN echo "Building app..." && cd app && LANDING_PAGE=/graph/extend PUBLIC_URL=https://playbook-workflow-builder.cloud npm run build

# TARGET: app -- production server with dependencies to run everything
FROM base AS app
USER root
RUN echo "Installing dev dependencies..." \
  && apt-get update -y \
  && apt-get install -y \
    libcurl4-gnutls-dev \
  && rm -rf /var/lib/apt/lists/*
USER ubuntu
RUN PYTHON_VERSION=3.11 R_VERSION=4.5.3 /install.sh
COPY --chown=ubuntu cli/setup.R app/setup.R
RUN echo "Running setup.R..." && R -e "source('app/setup.R')" && rm app/setup.R
COPY --from=prepare_requirements_txt_complete --chown=ubuntu /home/ubuntu/app app
RUN echo "Installing python dependencies..." && pip install -r app/requirements.txt && rm app/requirements.txt
COPY --from=prepare_build --chown=ubuntu /home/ubuntu/app app
WORKDIR /home/ubuntu/app
RUN chmod +x cli/wes-worker.sh cli/pwb.sh
ENV PORT=3000
ENV PORT=3005
ENV PYTHON_BIN=/home/ubuntu/.venv/bin/python
CMD ["npm", "start"]

# TARGET: app_minimal -- production server with dependencies to run just the webserver
FROM base AS app_minimal
COPY --from=prepare_build --chown=ubuntu /home/ubuntu/app app
WORKDIR /home/ubuntu/app
ENV PORT=3000
CMD ["npm", "run", "start"]

# TARGET: dev -- development environment with dependencies to run dev tools
FROM app AS dev
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
  && curl -LO https://github.com/ImageMagick/ImageMagick/archive/refs/tags/7.1.1-26.tar.gz \
  && tar xzf 7.1.1-26.tar.gz \
  && rm 7.1.1-26.tar.gz \
  && sh ./ImageMagick-7.1.1-26/configure --prefix=/usr/local --with-bzlib=yes --with-fontconfig=yes --with-freetype=yes --with-gslib=yes --with-gvc=yes --with-jpeg=yes --with-jp2=yes --with-png=yes --with-tiff=yes --with-xml=yes --with-gs-font-dir=yes --with-rsvg=yes \
  && make -j && make install && ldconfig /usr/local/lib/ \
  && apt-get clean \
  && apt-get autoremove \
  && rm -rf /var/lib/apt/lists/*
COPY --from=amacneil/dbmate /usr/local/bin/dbmate /usr/local/bin/dbmate
USER ubuntu
CMD ["/bin/bash"]
