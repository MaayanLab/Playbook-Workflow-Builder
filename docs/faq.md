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

## How can I add new icons to include in my components?

Icons are exported from the `icons/index.tsx`, two main types of icons are currently in use:
- icons from [Material Design Icons](https://materialdesignicons.com/), or from [BlueprintJS Icons](https://blueprintjs.com/docs/#icons). these are imported from their respective libraries and exported from `icons/mdi/index.ts` with a `path` (the svg path in the icon) and `title` a human readable label for that icon. By convention, they should end with `_icon`.
  ```ts
  // example material design icon for "gene"
  import { mdiDna } from '@mdi/js'; export const gene_icon = { path: mdiDna, title: 'Gene' }
  // example blueprintjs icon for "fork"
  export const fork_icon = { path: IconSvgPaths20.Fork.join(' '), title: 'Expand From this Step' }
  ```
- icons from some raster like a favicon. These icons can be added to `icons/services/src/*.png` and `npm run codegen:icons` can be used to trace those pngs into an svg path which gets inlined into `icons/services/index.ts`

Once registered, icons can be imported in your component like
```ts
import { some_icon } from '@/icons'
```

Note that the `Icon` type can be a list, the `<Icon>` react component in `app/components/icon.tsx` will display up to 4 different icons as one.

## Something in my development environment isn't working

We've Dockerized our dev environment, so if all else fails you should be able to use that given that you have `Docker` and `docker-compose` installed. It can be used with `docker-compose run dev`, in that shell you can execute any commands that weren't working.
