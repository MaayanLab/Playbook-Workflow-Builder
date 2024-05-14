import playbooks from '@/app/public/playbooksDemo'

console.log(
  [
    ['label', 'version', 'inputs', 'output', 'description'].join('\t'),
    ...playbooks.map(playbook => [
      playbook.label,
      playbook.version,
      playbook.inputs.map(input => input.meta.label).join('; '),
      playbook.outputs.map(output => output.meta.label).join('; '),
      JSON.stringify(playbook.description),
    ].join('\t')),
  ].join('\n')
)
