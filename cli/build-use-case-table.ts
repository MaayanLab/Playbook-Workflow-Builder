import playbooks from '@/app/public/playbooksDemo'

console.log(
  [
    ['label', 'version', 'inputs', 'output', 'description', 'url'].join('\t'),
    ...playbooks.map(playbook => [
      playbook.label,
      playbook.version,
      playbook.inputs.map(input => input.meta.label).join('; '),
      playbook.outputs.map(output => output.meta.label).join('; '),
      JSON.stringify(playbook.description),
      `https://playbook-workflow-builder.cloud/report/${playbook.id}`,
    ].join('\t')),
  ].join('\n')
)
