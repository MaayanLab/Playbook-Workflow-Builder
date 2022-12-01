# Frequently Asked Questions

## How do I build a codec for my data?

Besides expressing the shape of the data with [zod](https://zod.dev/), it's also possible to convert from existing schemas to zod. For example, if an OpenAPI/SmartAPI spec was already created, the schema of those parameters can be used with the [JSON Schema to Zod converter](https://stefanterdell.github.io/json-schema-to-zod-react/). Alternatively, if you just have raw JSON, you can get a head start using the [JSON to Zod converter](https://rsinohara.github.io/json-to-zod-react/).

## How can I safely pass private credentials for accessing an API?

The current way this can be done is through a service account token which can be passed as an environment variable to your component by including the token in the `.env` file and accessing it with `process.env.MY_SECRET_TOKEN`. Ultimately, this token will need to be shared with the current deployment maintainer to add to the deployment environment.

`.env`
```
#...
MY_SECRET_TOKEN=supersecrettoken
```

`components/[component]/index.ts`
```ts
MetaNode.createProcess()
  //...
  .resolve(async (props) => {
    const req = await fetch(`mypublic.api/route`, {
      headers: {
        // MY_SECRET_TOKEN is accessible here, only sever side code will have access to it
        Authorization: `Bearer ${process.env.MY_SECRET_TOKEN}`,
      },
    })
  })
```
