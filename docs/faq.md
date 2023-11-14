# Frequently Asked Questions

## How do I build a codec for my data?

Besides expressing the shape of the data with [zod](https://zod.dev/), it's also possible to convert from existing schemas to zod. For example, if an OpenAPI/SmartAPI spec was already created, the schema of those parameters can be used with the [JSON Schema to Zod converter](https://stefanterdell.github.io/json-schema-to-zod-react/). Alternatively, if you just have raw JSON, you can get a head start using the [JSON to Zod converter](https://rsinohara.github.io/json-to-zod-react/).

## How can I use files in MetaNodes?

To avoid storing large files in the database or in memory we pass around pointers to those files in the form of URLs. To interact with files you must first resolve it before you can manipulate it and potentially produce a new file which needs to be uploaded and turned back into a pointer. Since how this is done is subject to change (i.e. storing on the filesystem, in s3, or other mechanisms), helper functions are available for managing files.

### Access files from Python-Implemented Components

```python
# turning a file object into content:
#  when reading sequentially, you can read it as a stream directly from
#  wherever it is. it behaves like a file object
from components.core.file import file_as_stream
with file_as_stream(file) as fh: # you get something like a file handle here
  for line in fh:
    print(line)

# if your file needs seek support (like some database formats, i.e. h5)
#  you'll need to use file_as_path which ensures its on the disk
from components.core.file import file_as_path
with file_as_path(file['url']) as path: # you get an actual path rather than a handle
  with open(path, 'r+') as fh:
    print(fh.read())

# turning content into a file object
from components.core.file import upsert_file
with upsert_file('.ext') as f:
  # f.file is a python file handler
  f.file.write('test')
# after the context manager, f.url contains the uploaded file url
return f
```

### Access files from NodeJS-Implemented Resolvers

Please note only **resolvers** should be operating on files, *views* and *prompts* should operate on prepared datatypes which resolvers could set up. That-is picking options in a way that depends on information in a file would need to prepare those options in a resolver. This is because potentially large files should not be downloaded to an end user's browser which is where views and prompts run.

```tsx
import { fileAsStream } from  '@/components/core/file/api/download'
import { fileFromStream } from  '@/components/core/file/api/upload'
import FormData from 'form-data'
import axios from 'axios'

export const SomeFileOp = MetaNode('SomeFileOp')
  .meta({ label: 'Some File Operation', description: 'Perform an operation' })
  .inputs({ file: FileURL })
  .output(FileURL)
  .resolve(async (props) => {
    // read the file as a stream
    const fileReader: any = await fileAsStream(props.inputs.file)
    // send that file somewhere for processing
    const formData = new FormData()
    formData.append('file', fileReader, props.inputs.file.filename)
    const res = await axios.post('https://example.com/upload', formData, {
      headers: { ...formData.getHeaders() },
      responseType: 'stream',
    })
    // create a new file object from a stream
    const file = await fileFromStream(res.data, `derived.${props.inputs.file.filename}`)
    return file
  })
  .story(props => `The performed some operation on the file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''}.`)
  .build()
```

## How can I safely pass private credentials for accessing an API?

The current way this can be done is through a service account token which can be passed as an environment variable to your component by including the token in the `.env` file and accessing it with `process.env.MY_SECRET_TOKEN`. Ultimately, this token will need to be shared with the current deployment maintainer to add to the deployment environment.

`.env`
```
#...
MY_SECRET_TOKEN=supersecrettoken
```

`components/[component]/index.ts`
```ts
MetaNode('MyProcess')
  //...
  .resolve(async (props) => {
    const req = await fetch(`mypublic.api/route`, {
      headers: {
        // MY_SECRET_TOKEN is accessible here, only sever side code will have access to it
        Authorization: `Bearer ${process.env.MY_SECRET_TOKEN}`,
      },
    })
  })
  //...
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

## Learn More

[Find other topics in the Playbook Workflow Builder Developer Guide](./index.md).
