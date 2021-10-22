import * as React from "react";
import { GenericPlotSelectOption, IPlotSelectExtension } from "./PlotSelect";

export const PlotSelectJaccardExtension: IPlotSelectExtension = {
  component: ({
    option,
    setOption,
    availableNearestNeighbors,
  }: {
    availableNearestNeighbors: string[];
  } & GenericPlotSelectOption<string>) => {
    const [key, jaccardEmb1, jaccardEmb2, jaccardNN] = option?.split("=") || [];
    const isJaccard = key === "jaccard";

    return isJaccard ? (
      <>
        <div className="col">
          <select
            className="form-control"
            value={jaccardEmb1 || ""}
            onChange={(e) =>
              setOption(`jaccard=${e.currentTarget.value}=${jaccardEmb2 || ""}=${jaccardNN || ""}` as any)
            }
          >
            <option value="">Select...</option>
            {availableNearestNeighbors.map((nn) => (
              <option value={nn}>{nn}</option>
            ))}
          </select>
        </div>
        <div className="col">
          <select
            className="form-control"
            value={jaccardEmb2 || ""}
            onChange={(e) =>
              setOption(`jaccard=${jaccardEmb1 || ""}=${e.currentTarget.value}=${jaccardNN || ""}` as any)
            }
          >
            <option value="">Select...</option>
            {availableNearestNeighbors.map((nn) => (
              <option value={nn}>{nn}</option>
            ))}
          </select>
        </div>
        <div className="col">
          <input
            className="form-control"
            type="number"
            min={1}
            max={50}
            value={+jaccardNN || ""}
            onChange={(e) =>
              setOption(`jaccard=${jaccardEmb1 || ""}=${jaccardEmb2 || ""}=${e.currentTarget.value}` as any)
            }
          />
        </div>
      </>
    ) : null;
  },
  additionalOptions: [
    {
      label: "Ratio of matching nearest neighbors",
      value: "jaccard=",
    },
  ],
};
