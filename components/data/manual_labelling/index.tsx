import React from 'react';
import { MetaNode } from '@/spec/metanode';
import { z } from 'zod';

export const string1 = MetaNode(`string1`)
  .meta({
    label: 'Simple String1 MetaNode',
    description: 'A MetaNode that represents a simple string1 value',
  })
  .codec(z.string())
  .view(term => {
    // Replace newline characters with HTML line breaks
    const formattedTerm = term.replace(/\n/g, '<br>');

    return <div dangerouslySetInnerHTML={{ __html: formattedTerm }} />;
  })
  .build();

export const ManualLabelling = MetaNode('ManualLabelling')
  .meta({
    label: 'Manual Labelling of Samples',
    description: 'Manually label samples as either control or perturbation.',
  })
  .inputs()
  .output(string1)
  .prompt(props => {
    const [tableData, setTableData] = React.useState([
      ['', ''],
      ['', ''],
      ['', ''],
      ['', ''],
      ['', ''],
      ['', ''],
      ['', ''],
      ['', ''],
      ['', ''],
      ['', ''],
    ]);
    //just set to 10 rows for now because this is based on how many samples exist 

    const handleCellChange = (rowIndex: number, columnIndex: number, value: string) => {
      const newTableData = [...tableData];
      newTableData[rowIndex][columnIndex] = value;
      setTableData(newTableData);
    };

    const handleCellSubmit = () => {
      const controlSamples: string[] = [];
      const perturbationSamples: string[] = [];

      // Separate control and perturbation samples
      tableData.forEach(row => {
          if (row[1] === 'Control') {
            controlSamples.push(row[0]);
          } else if (row[1] === 'Perturbation') {
            perturbationSamples.push(row[0]);
        }
      });

      // Create the final output string with line breaks
      const outputString =
        `{"controls": [${controlSamples.map(sample => `"${sample}"`).join(', ')}],\n` +
        ` "perturbations": [${perturbationSamples.map(sample => `"${sample}"`).join(', ')}]}`;

      props.submit(outputString);
    };

    return (
      <div>
        <table style={{ borderCollapse: 'collapse', border: '1px solid black' }}>
          <thead>
            <tr>
              <th style={{ border: '1px solid black', padding: '5px' }}>Samples</th>
              <th style={{ border: '1px solid black', padding: '5px' }}>Type: Control or Perturbation</th>
            </tr>
          </thead>
          <tbody>
            {tableData.map((row, rowIndex) => (
              <tr key={rowIndex}>
                {row.map((cellValue, columnIndex) => (
                  <td
                    key={columnIndex}
                    style={{
                      border: '1px solid black',
                      padding: '5px',
                    }}
                  >
                    {columnIndex === 1 ? (
                      <select
                        value={cellValue}
                        onChange={evt =>
                          handleCellChange(rowIndex, columnIndex, evt.target.value)
                        }
                      >
                        <option value="">Select</option>
                        <option value="Control">Control</option>
                        <option value="Perturbation">Perturbation</option>
                      </select>
                    ) : (
                      <input
                        type="text"
                        value={cellValue}
                        onChange={evt =>
                          handleCellChange(rowIndex, columnIndex, evt.target.value)
                        }
                        style={{
                          width: '100%',
                          boxSizing: 'border-box',
                          border: 'none',
                          outline: 'none',
                          textAlign: 'center',
                        }}
                      />
                    )}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
        <button onClick={handleCellSubmit}>Submit</button>
      </div>
    );
  })
  .story(props =>
    'The samples were then labelled as either control or perturbation to allow for further analysis.'
  )
  .build();




