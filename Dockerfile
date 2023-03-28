FROM node:latest as prepare_package_json
WORKDIR /app
COPY . /app
RUN find /app -type d -name "node_modules" -exec rm -rf {} + \
  && find /app -type f -a \! \( -name "package.json" -o -name "package-lock.json" \) -delete \
  && find /app -type d -empty -delete

FROM node:latest as prepare_requirements_txt
WORKDIR /app
COPY . /app
RUN find /app -type d -name "node_modules" -exec rm -rf {} + \
  && find /app -type f -a \! \( -name "requirements.txt" -o -name "enumerate-requirements.ts" \) -delete \
  && find /app -type d -empty -delete

FROM node:latest as prepare_npm_i
WORKDIR /app
COPY --from=prepare_package_json /app /app
RUN echo "Installing NodeJS dependencies..." && npm i

FROM node:latest as prepare_requirements_txt_complete
WORKDIR /app
COPY --from=prepare_npm_i /app /app
COPY --from=prepare_requirements_txt /app /app
RUN npm run codegen:requirements \
  && mv /app/requirements.txt /tmp/requirements.txt \
  && rm -r /app \
  && mkdir /app \
  && mv /tmp/requirements.txt /app

FROM node:latest as prepare_build
WORKDIR /app
COPY --from=prepare_npm_i /app /app
COPY . /app
RUN echo "Building app..." && LANDING_PAGE=/graph/extend npm run build

FROM node:latest as app
WORKDIR /app
ENV PYTHON_BIN="python3"
RUN echo "Installing python..." && apt-get -y update && apt-get -y install python3-dev python3-pip && rm -rf /var/lib/apt/lists/*
COPY --from=prepare_requirements_txt_complete /app /app
RUN echo "Installing python dependencies..." && pip install -r /app/requirements.txt && rm /app/requirements.txt
COPY --from=prepare_build /app /app
EXPOSE 3000
CMD ["npm", "run", "start"]
