import * as React from "react";
import { GenericPlotSelectOption, IPlotSelectExtension } from "./PlotSelect";

export const PlotSelectEvalExtension: IPlotSelectExtension = {
  component: ({ option, setOption }: {} & GenericPlotSelectOption<string>) => {
    const [key, func] = option?.split("=") || [];
    const isEval = key === "eval";

    return isEval ? (
      <>
        <div className="col-8">
          <div className="input-group">
            <span className="input-group-text">
              <code>(p) ={">"} </code>
            </span>
            <input
              className="form-control"
              type="text"
              value={func || ""}
              onChange={(e) => setOption(`eval=${e.currentTarget.value || ""}` as any)}
            />
          </div>
        </div>
      </>
    ) : null;
  },
  additionalOptions: [
    {
      label: "Custom function",
      value: "eval=",
    },
  ],
};
