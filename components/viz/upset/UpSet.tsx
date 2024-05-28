/**
 * Credit to Stephanie Olaiya -- from G2SG
 */
import * as d3 from "d3";
import React from "react"
import { UpsetTooltip, UpsetInteractionData } from "./UpSetTooltip";

type UpsetData = {
    name: string,
    values: string[]
}[]

type SoloIntersectionType = {
    name: string,
    setName: string,
    num: number
    values: string[]
}


// Upset plot code adapted from https://github.com/chuntul/d3-upset
const formatIntersectionData = (data: UpsetData) => {
    // compiling solo set data - how many values per set
    const soloSets: SoloIntersectionType[] = [];

    // nameStr is for the setName, which makes it easy to compile
    // each name would be A, then B, so on..
    const nameStr = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.substring(0, data.length);
    data.forEach((x, i) => {
        soloSets.push({
            name: x.name,
            setName: nameStr.substring(i, i+1),
            num: x.values.length,
            values: x.values,
        });
    });

    // compiling list of intersection names recursively
    // ["A", "AB", "ABC", ...]
    const getIntNames = (start: number, end: number, nameStr: string): string[] => {
        // eg. BCD
        const name = nameStr.substring(start, end);
        // when reaching the last letter
        if (name.length === 1) {
            return [name];
        }
        const retArr = getIntNames(start + 1, end, nameStr);
        // eg. for name = BCD, would return [B] + [BC,BCD,BD] + [C,CD,D]
        return [name[0]].concat(retArr.map((x) => name[0] + x), retArr);
    };

    let intNames = getIntNames(0, nameStr.length, nameStr);

    // removing solo names
    intNames = intNames.filter((x) => x.length !== 1);

    let intersections: SoloIntersectionType[] = [];

    // compile intersections of values for each intersection name
    intNames.forEach((intName) => {
        // collecting all values: [pub1arr, pub2arr, ...]
        const values = intName.split('').map((set) => {
            const sets = soloSets.find((x) => x.setName === set)?.values
            if (sets) {
                return sets
            } else { return [] }
        }
        );

        // getting intersection
        // https://stackoverflow.com/questions/37320296/how-to-calculate-intersection-of-multiple-arrays-in-javascript-and-what-does-e
        const result = values.reduce((a, b) => a.filter((c) => b.includes(c)));
        intersections.push({
            name: intName.split('').map((set) => soloSets.find((x) => x.setName === set)?.name).join(' + '),
            setName: intName,
            num: result.length,
            values: result,
        });
    });

    // taking out all 0s
    intersections = intersections.filter((x) => x.num !== 0); // changed from .value 
    return { intersections, soloSets };
};

// include solo sets with all its data
const insertSoloDataAll = (intersections: SoloIntersectionType[], soloSets: SoloIntersectionType[]) => {
    soloSets.forEach(x => {
        intersections.push(x);
    });
    return intersections;
};

// include solo sets with only the values that ARE NOT in other sets
const insertSoloDataOutersect = (intersections: SoloIntersectionType[], soloSets: SoloIntersectionType[]) => {
    soloSets.forEach(x => {
        // compile all unique values from other sets except current set
        const otherSets = Array.from(new Set(soloSets.map(y => y.setName === x.setName ? [] : y.values).flat()));

        // subtract otherSets values from current set values
        const values = x.values.filter(y => !otherSets.includes(y));
        intersections.push({
            name: x.name,
            setName: x.setName,
            num: values.length,
            values: values,
        })

    })
    return intersections;
}




export function UpsetPlotV2({ selectedSets, setOverlap }: {
    selectedSets: ({
        alphabet: string;
        genes: { gene_symbol: string }[];
    } & { alphabet: string })[] | undefined;
    setOverlap: React.Dispatch<React.SetStateAction<{name: string, overlapGenes: string[]}>>;

}) {
    const [hoveredCell, setHoveredCell] = React.useState<UpsetInteractionData | null>(null);

    const { data, soloSets }: { data: SoloIntersectionType[]; soloSets: SoloIntersectionType[]; } = React.useMemo(() => {
        if (!selectedSets) {
            const data: SoloIntersectionType[] = []
            const soloSets: SoloIntersectionType[] = []
            return { data, soloSets };
        }

        const setData = selectedSets.map((geneset, i) => { return { name: geneset.alphabet, values: geneset.genes.map((gene) => gene.gene_symbol) } })

        // calculating intersections WITHOUT solo sets
        const { intersections, soloSets } = formatIntersectionData(setData);

        // putting the solo sets in:
        // to include solo sets with all its data, call this function
        const data = insertSoloDataAll(intersections, soloSets);

        // to include solo sets with only the values that ARE NOT in other sets
        // ie. the outersect of values in the solo sets, call this function
        // export const allData = insertSoloDataOutersect(intersections, soloSets);
        return { data, soloSets }
    }, [selectedSets])

    const allSetNames = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.substring(0, soloSets.length).split('');

    // position and dimensions
    const margin = {
        top: 30,
        right: 0,
        bottom: Math.max(250, 35 * soloSets.length),
        left: 50,
    };
    const width = (30 * data.length) + 100 + 60 + 200;
    // console.log(soloSets.length)
    const height = Math.max(400, 60 * soloSets.length);

    // The bounds (=area inside the axis) is calculated by substracting the margins
    const boundsWidth = width - margin.right - margin.left;
    const boundsHeight = height - margin.top - margin.bottom;

    // sort data decreasing
    data.sort((a, b) => b.num - a.num);

    const nums = data.map((x) => x.num);

    // set range for data by domain, and scale by range
    const xrange = d3.scaleLinear().range([0, boundsWidth]).domain([0, data.length])

    const yrange = d3
        .scaleLinear()
        .range([boundsHeight, 0])
        .domain([0, d3.max(nums) as number])

    // left yaxis
    const range = yrange.range();
    const pixelsPerTick = 20 // increase spacing between ticks
    const lineheight = range[0] - range[1];
    const numberOfTicksTarget = Math.floor(lineheight / pixelsPerTick);
    const ticks = yrange.ticks(numberOfTicksTarget).map((value) => ({
        value,
        yOffset: yrange(value),
    }));
    const leftAxis =
        <>
            {/* Main vertical line */}
            <path
                d={["M", 0, range[0], "L", 0, range[1]].join(" ")}
                fill="none"
                stroke="currentColor"
            />
            {/* Axis label*/}
            <text
                x={(margin.bottom - height) / 2}
                y={range[1] - margin.left}
                textAnchor="middle"
                dominantBaseline="middle"
                fill='black'
                fontSize={12}
                transform="rotate(-90)"
            >
                Intersection Size
            </text>
            {/* Ticks and labels */}
            {ticks.map(({ value, yOffset }) => (
                <g key={value} transform={`translate(0, ${yOffset})`}>
                    <line x2={-6} stroke="currentColor" />
                    <text
                        key={value.toString() + '-1'}
                        style={{
                            fontSize: "10px",
                            textAnchor: "middle",
                            transform: "translateX(-20px)",
                        }}
                    >
                        {value}
                    </text>
                </g>
            ))}
        </>

    const rad = 10;
    const bottomAxisrange = xrange.range();
    const bottomAxis =
        <>
            {/* Main horizontal line */}
            <path
                d={["M", bottomAxisrange[0], range[0], "L", 9 + (data.length) * (rad * 2.7), range[0]].join(" ")}
                fill="none"
                stroke="currentColor"
            />
        </>



    const labels = soloSets.map((set, i) => {
        return (
            <text
                key={i}
                x={-30}
                y={5 + i * (rad * 2.7)}
                textAnchor="end"
                dominantBaseline="middle"
                fill='black'
                fontSize={10}
                onMouseEnter={(e) => {
                    setHoveredCell({
                        setLabel: set.name,
                        xPos: 9 + i * (rad * 2.7) + margin.left,
                        yPos: boundsHeight + yrange(set.num) + margin.top,
                        value: set.num,
                    });
                }}
                onMouseLeave={() => setHoveredCell(null)}
                cursor="pointer"
            >
                {set.setName}
            </text>
        );
    });

    const soloSetsLengthNums = soloSets.map((x) => x.num);
    // set range for data by domain, and scale by range
    const xrangeSet = d3.scaleLinear().range([0, 100]).domain([0, d3.max(soloSetsLengthNums) as number])
    const invertedXrangeSet = d3.scaleLinear()
        .range([100, 0]) // Reverse the range
        .domain(xrangeSet.domain()); // Use the same domain as the original scale 
    const setBars = soloSets.map((d, i) => {
        let fillColor
        if (d.name == hoveredCell?.setLabel) {
            fillColor = '#FFC000';

        } else {
            // fillColor = '#02577b';
            fillColor = '#000000';
        }
        return (
            <rect
                key={i}
                r={4}
                // x={0}
                y={i * (rad * 2.7)}
                x={invertedXrangeSet(d.num) + margin.left}
                width={Math.abs(xrangeSet(d.num))}
                height={(rad * 2.7) - 9}
                opacity={1}
                fill={fillColor}
                stroke={"white"}
                onMouseEnter={(e) => {
                    setHoveredCell({
                        setLabel: d.name,
                        xPos: invertedXrangeSet(d.num) + 60 + margin.left,
                        yPos: boundsHeight + i * (rad * 2.7) + margin.top + 30,
                        value: d.num,
                    });
                }}
                onMouseLeave={() => setHoveredCell(null)}
                cursor="pointer"
                onMouseDown={() => setOverlap({name: d.setName, overlapGenes: d.values})}
            />

        );
    });

    const setSizerange = invertedXrangeSet.range();
    const setSizePixelsPerTick = 20 // increase spacing between ticks
    const setSizeLineheight = setSizerange[0] - setSizerange[1];
    const setSizeNumberOfTicksTarget = Math.floor(setSizeLineheight / setSizePixelsPerTick);
    const setSizeTicks = invertedXrangeSet.ticks(setSizeNumberOfTicksTarget).map((value) => ({
        value,
        xOffset: invertedXrangeSet(value),
    }));
    const soloBottomAxisrange = invertedXrangeSet.range();
    const soloBottomAxis =
        <>
            <path
                d={["M", soloBottomAxisrange[0] + margin.left, soloSets.length * (rad * 2.7), "L", soloBottomAxisrange[1] + margin.left, soloSets.length * (rad * 2.7)].join(" ")}
                fill="none"
                stroke="currentColor"
            />
            <text
                x={soloBottomAxisrange[0] / 2 + margin.left}
                y={soloSets.length * (rad * 2.7) + 30}
                textAnchor="middle"
                dominantBaseline="middle"
                fill='black'
                fontSize={12}

            >
                Set Size
            </text>
            {/* Ticks and labels */}
            {setSizeTicks.map(({ value, xOffset }) => (
                <g key={value} transform={`translate(${xOffset + margin.left}, ${soloSets.length * (rad * 2.7) + 6})`}>
                    <line y2={-6} stroke="currentColor" />
                    <text
                        key={value.toString() + '-2'}
                        style={{
                            fontSize: "10px",
                            textAnchor: "middle",
                            transform: `translate(0px, 10px) rotate(45deg)`,
                        }}
                    >
                        {value}
                    </text>
                </g>
            ))}
        </>

    // bars 
    const bars = data.map((d, i) => {
        const x = xrange(9 + i * (rad * 2.7));
        const y = yrange(yrange(d.num));

        if (d.values === null || !x || !y) {
            return;
        }
        let fillColor
        if (d.name == hoveredCell?.setLabel) {
            fillColor = '#FFC000';

        } else {
            fillColor = '#000000';
        }
        return (
            <rect
                key={i}
                r={4}
                x={9 + i * (rad * 2.7)}
                y={yrange(d.num)}
                width={20}
                height={boundsHeight - yrange(d.num)}
                opacity={1}
                fill={fillColor}
                // rx={5}
                stroke={"white"}
                onMouseEnter={(e) => {
                    setHoveredCell({
                        setLabel: d.name,
                        xPos: 9 + i * (rad * 2.7) + margin.left + 200,
                        yPos: yrange(d.num) + margin.top,
                        value: d.num,
                    });

                }}
                onMouseLeave={() => setHoveredCell(null)}
                cursor="pointer"
                onMouseDown={() => setOverlap({name: d.setName, overlapGenes: d.values})}
            />
        );
    })

    // circles for image
    const circles = data.map((x, i) => {
        return allSetNames.map((y, j) => {
            let fillColor
            if (x.setName.indexOf(y) !== -1) {
                fillColor = '#000000';
            } else {
                fillColor = 'silver';
            }

            return (
                <circle
                    key={i + ',' +  j}
                    r={rad}
                    cx={i * (rad * 2.7)}
                    cy={j * (rad * 2.7)}
                    opacity={1}
                    fill={fillColor}
                />
            );
        })
    })

    // add lines to circles
    const lines = data.map((x, i) => {
        return allSetNames.map((y, j) => {
            return (
                <line
                    key={i +',' +j + '-circle'}
                    id={`setline${i}`}
                    x1={i * (rad * 2.7)}
                    x2={i * (rad * 2.7)}
                    y1={allSetNames.indexOf(x.setName[0]) * (rad * 2.7)}
                    y2={allSetNames.indexOf(x.setName[x.setName.length - 1]) * (rad * 2.7)}
                    stroke="#000000"
                    strokeWidth={4}
                />
            );
        })
    })


    return (
        <div style={{ position: "relative", overflow: "auto", maxWidth: 700 }} >
            <svg width={width} height={height} id='svg'>
                <g
                    width={boundsWidth}
                    height={boundsHeight}
                    transform={`translate(${[margin.left, 0].join(",")})`}
                >
                    <g id='upsetBars'
                        transform={`translate(${200},${margin.top})`}
                    >
                        <g id='chart'
                            transform={'translate(1,0)'}
                        >
                            {bars}

                        </g>
                        {leftAxis}
                        {bottomAxis}

                    </g>

                    <g id='setBars'
                        transform={`translate(${[0, boundsHeight + margin.top + 25].join(",")})`}
                    >
                        {setBars}
                        {soloBottomAxis}
                    </g>

                    <g id='upsetCircles'
                        transform={`translate(${[220, boundsHeight + margin.top + 30].join(",")})`} // change 100 to another value
                    >
                        {labels}
                        {circles}
                        {lines}
                    </g>
                </g>
            </svg>
            <UpsetTooltip interactionData={hoveredCell} width={width} height={height} />
        </div>
    );

}