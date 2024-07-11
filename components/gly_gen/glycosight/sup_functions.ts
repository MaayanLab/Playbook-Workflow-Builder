import { GlycoSightOutputType, GlycoSightFileURLType } from "./data_models";
import { GlycoSightOutputNode } from ".";
import { z } from "zod";

const devGSURL = process.env.NODE_ENV === "development" ? "http://localhost:5000/api/" : "https://aws.glygen.org/glycosight/api/"

export async function TestGlycoSightAPI(url: string) : Promise<GlycoSightOutputType> {
    // const res = await fetch(`${devGSURL}/dummy-upload?q=${url}&n=${fileName}`, { method: "POST" })
    const res = await fetch(`${devGSURL}test-access?q=${url}`, { method: "POST" });
    const dummyResult = res.json();
    return dummyResult;
}

export async function UploadAndAnalyze(file: GlycoSightFileURLType, fileData: Buffer) : Promise<GlycoSightOutputType> {
    const res = await fetch(`${devGSURL}upload-and-analyze?n=${file.filename}`, 
        { 
            method: "POST", 
            body: fileData, 
        //     headers: {
        //         "Content-length": fileData.length
        // } 
        }
    );
    const result = await res.json();

    // const dummy = {
    //     results: [
    //         {
    //             UniProtAcc: "12345",
    //             AAPosition: 42,
    //             SpectralCount: 1,
    //             DistinctPeptideCount: 1,
    //             Peptides: "foobar",
    //         },
    //     ]
    // }
    
    return result;
}


