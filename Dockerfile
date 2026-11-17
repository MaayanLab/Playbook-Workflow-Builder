# NodeJS, from its official image: copied into the base (below) rather than installed with nvm,
#  which extracts it with tar: that fails when building for amd64 on Apple Silicon (emulated, tar's
#  "Cannot open: Function not implemented"), copying doesn't run in the emulated container
FROM node:26-slim AS node

FROM maayanlab/base:1.6.0 AS base
USER root
ENV PUPPETEER_SKIP_DOWNLOAD=true
RUN echo "Installing system dependencies (git+puppeteer deps)..." \
  && apt-get update \
  && apt-get -y install \
    biber \
    build-essential \
    curl \
    git \
    gnupg \
    latexmk \
    libatomic1 \
    python3 \
    texlive-bibtex-extra \
    texlive-fonts-extra \
    texlive-latex-extra \
    texlive-luatex \
  && install -d -m 0755 /etc/apt/keyrings \
  && curl --location --silent --fail https://dl-ssl.google.com/linux/linux_signing_key.pub | gpg --dearmor -o /etc/apt/keyrings/google-chrome.gpg \
  && echo "deb [arch=amd64 signed-by=/etc/apt/keyrings/google-chrome.gpg] http://dl.google.com/linux/chrome/deb/ stable main" > /etc/apt/sources.list.d/google-chrome.list \
  && apt-get update \
  && apt-get install google-chrome-stable -y --no-install-recommends \
  && rm -rf /var/lib/apt/lists/*
ENV PUPPETEER_EXECUTABLE_PATH="/usr/bin/google-chrome"
COPY --from=node /usr/local/bin/ /usr/local/bin/
COPY --from=node /usr/local/lib/node_modules/ /usr/local/lib/node_modules/
COPY --from=node /usr/local/include/node/ /usr/local/include/node/
USER ubuntu

# The JavaScript is prepared & built on the build machine's own platform ($BUILDPLATFORM): its output is
#  the same on every platform, and building it emulated (e.g. amd64 on Apple Silicon) is many times
#  slower. Only the runtime's node_modules (which have native binaries) are installed for the target.
FROM --platform=$BUILDPLATFORM node:26-slim AS build_base
# (to compile native dependencies, e.g. better-sqlite3)
RUN apt-get update && apt-get -y install --no-install-recommends python3 make g++ && rm -rf /var/lib/apt/lists/*
ENV PUPPETEER_SKIP_DOWNLOAD=true
WORKDIR /home/ubuntu

FROM build_base AS prepare_src
COPY . app

# (only the package.json files: the dependencies are installed again only when they change)
FROM prepare_src AS prepare_package_json
RUN find app -type d -name "node_modules" -exec rm -rf {} + \
  && find app -type f -a \! \( -name "package.json" -o -name "package-lock.json" -o -name ".puppeteerrc.cjs" \) -delete \
  && find app -type d -empty -delete

FROM prepare_src AS prepare_requirements_txt
RUN find app -type d -name "node_modules" -exec rm -rf {} + \
  && find app -type f -a \! \( -name "requirements.txt" -o -name "enumerate-requirements.ts" \) -delete \
  && find app -type d -empty -delete

# the dependencies for building (on the build platform), npm's downloads are kept between builds
FROM build_base AS prepare_npm_i_build
COPY --from=prepare_package_json /home/ubuntu/app app
RUN --mount=type=cache,target=/root/.npm \
  echo "Installing NodeJS dependencies (build)..." && cd app && npm i

FROM prepare_npm_i_build AS prepare_requirements_txt_complete
COPY --from=prepare_requirements_txt /home/ubuntu/app app
RUN cd app && npm run codegen:requirements

# Next's cache is kept between builds (& out of the image)
FROM prepare_npm_i_build AS prepare_build
COPY --from=prepare_src /home/ubuntu/app app
RUN --mount=type=cache,target=/home/ubuntu/app/app/.next/cache \
  echo "Building app..." && cd app && npm run build

# the dependencies for running (on the target platform)
FROM base AS prepare_npm_i
COPY --from=prepare_package_json --chown=ubuntu /home/ubuntu/app app
RUN --mount=type=cache,target=/home/ubuntu/.npm,uid=1000,gid=1000 \
  echo "Installing NodeJS dependencies..." && cd app && npm i

# TARGET: app -- production server with dependencies to run everything
#  (ordered from what changes least to most: python, node_modules, the source, the build)
FROM base AS app
USER ubuntu
RUN PYTHON_VERSION=3.13 /install.sh
COPY --from=prepare_requirements_txt_complete --chown=ubuntu /home/ubuntu/app/requirements.txt app/requirements.txt
RUN --mount=type=cache,target=/home/ubuntu/.cache/uv,uid=1000,gid=1000 \
  echo "Installing python dependencies..." && UV_LINK_MODE=copy uv pip install -r app/requirements.txt && rm app/requirements.txt
COPY --from=prepare_npm_i --chown=ubuntu /home/ubuntu/app app
COPY --from=prepare_src --chown=ubuntu /home/ubuntu/app app
COPY --from=prepare_build --chown=ubuntu /home/ubuntu/app/app/.next app/app/.next
WORKDIR /home/ubuntu/app
RUN chmod +x cli/wes-worker.sh cli/pwb.sh
ENV PORT=3000
ENV PYTHON_BIN=/home/ubuntu/.venv/bin/python
CMD ["npm", "start"]

# TARGET: app_minimal -- production server with dependencies to run just the webserver
FROM base AS app_minimal
COPY --from=prepare_npm_i --chown=ubuntu /home/ubuntu/app app
COPY --from=prepare_src --chown=ubuntu /home/ubuntu/app app
COPY --from=prepare_build --chown=ubuntu /home/ubuntu/app/app/.next app/app/.next
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
USER ubuntu
CMD ["/bin/bash"]
