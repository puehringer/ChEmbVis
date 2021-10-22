import * as React from "react";
import { GenericPlotSelectOption, IPlotSelectExtension } from "./PlotSelect";

export const PlotSelectClusterExtension: IPlotSelectExtension = {
  component: ({
    option,
    setOption,
    availableNearestNeighbors,
    availableClusters,
    availableProperties,
  }: {
    availableNearestNeighbors: string[];
    availableClusters: string[];
    availableProperties: string[];
  } & GenericPlotSelectOption<string>) => {
    const [key, emb] = option?.split("=") || [];
    const isNNDiff = key === "cluster";

    return isNNDiff ? (
      <>
        <div className="col">
          <select
            className="form-control"
            value={emb || ""}
            onChange={(e) =>
              setOption(`cluster=${e.currentTarget.value}` as any)
            }
          >
            <option value="">Select...</option>
            {availableClusters.map((nn) => (
              <option value={nn}>{nn}</option>
            ))}
          </select>
        </div>
      </>
    ) : null;
  },
  additionalOptions: [
    {
      label: "Clusters",
      value: "cluster=",
    },
  ],
};
