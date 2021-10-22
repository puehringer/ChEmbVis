import * as React from "react";
import { IEnabledProjection } from "../../interfaces";
import { PlotSelect } from "../PlotSelect";
import uniqueId from "lodash.uniqueid";
import { PlotSelectJaccardExtension } from "../PlotSelectJaccardExtension";
import { PlotSelectNNDiffExtension } from "../PlotSelectNNDiffExtension";
import { PlotSelectKNNExtension } from "../PlotSelectKNNExtension";
import { PlotSelectClusterExtension } from "../PlotSelectClusterExtension";
import { PlotSelectEvalExtension } from "../PlotSelectEvalExtension";

export const ProjectionSettingsForm = ({
  inline,
  config,
  setConfig,
  availableProjections,
  availableNearestNeighbors,
  availableClusters,
  availableProperties,
  availableOpacityProperties,
  availableConnectByProperties,
}: {
  inline: boolean;
  config: IEnabledProjection;
  setConfig(value: IEnabledProjection): void;
  availableProjections?: string[];
  availableNearestNeighbors: string[];
  availableClusters: string[];
  availableProperties: string[];
  availableOpacityProperties: string[];
  availableConnectByProperties: string[];
}) => {
  const [id, setId] = React.useState<string>(uniqueId());
  return (
    <form className="row">
      {availableProjections ? (
        <PlotSelect
          multi={false}
          inline={inline}
          id={`projection${id}Select`}
          label="Projection"
          option={config.projection}
          options={availableProjections}
          setOption={(projection: string) => {
            setConfig({ ...config, projection });
          }}
          availableNearestNeighbors={availableNearestNeighbors}
          availableClusters={availableClusters}
          availableProperties={availableProperties}
        />
      ) : null}
      <PlotSelect
        multi={false}
        inline={inline}
        id={`colorEncoding${id}Select`}
        label="Color Encoding"
        // disabled={groupBy !== null}
        option={config.plotOptions.colorBy}
        options={availableProperties}
        setOption={(colorCoding: string) => {
          setConfig({ ...config, plotOptions: { ...config.plotOptions, colorBy: colorCoding } });
        }}
        availableNearestNeighbors={availableNearestNeighbors}
        availableClusters={availableClusters}
        availableProperties={availableProperties}
        extensions={[PlotSelectJaccardExtension, PlotSelectNNDiffExtension, PlotSelectClusterExtension, PlotSelectEvalExtension]}
      />
      <PlotSelect
        multi={false}
        inline={inline}
        id={`opacity${id}Select`}
        label="Opacity By"
        // disabled={colorCoding !== null}
        option={typeof config.plotOptions.opacityBy === "number" ? "constant" : config.plotOptions.opacityBy}
        options={availableOpacityProperties}
        setOption={(opacityBy: string) => {
          setConfig({
            ...config,
            plotOptions: { ...config.plotOptions, opacityBy: opacityBy === "constant" ? 0.5 : opacityBy },
          });
        }}
        availableNearestNeighbors={availableNearestNeighbors}
        availableClusters={availableClusters}
        availableProperties={availableProperties}
        extensions={[PlotSelectJaccardExtension, PlotSelectNNDiffExtension, PlotSelectClusterExtension, PlotSelectEvalExtension]}
      />
      {typeof config.plotOptions.opacityBy === "number" ? (
        <div className="mb-3 me-sm-2">
          {/* <label for="formControlRange" className="form-label">Example Range input</label> */}
          <input
            type="range"
            className="form-range"
            min={0}
            max={1}
            step={0.01}
            value={config.plotOptions.opacityBy}
            onChange={(e) => {
              setConfig({
                ...config,
                plotOptions: { ...config.plotOptions, opacityBy: e.currentTarget.valueAsNumber },
              });
            }}
          />
        </div>
      ) : null}
      <PlotSelect
        multi={false}
        inline={inline}
        id={`grouping${id}Select`}
        label="Group By"
        // disabled={colorCoding !== null}
        option={config.plotOptions.groupBy}
        options={availableProperties}
        setOption={(groupBy: string) => {
          setConfig({ ...config, plotOptions: { ...config.plotOptions, groupBy } });
        }}
        availableNearestNeighbors={availableNearestNeighbors}
        availableClusters={availableClusters}
        availableProperties={availableProperties}
        extensions={[PlotSelectJaccardExtension, PlotSelectNNDiffExtension, PlotSelectClusterExtension, PlotSelectEvalExtension]}
      />
      <PlotSelect
        multi={false}
        inline={inline}
        id={`sizeBy${id}Select`}
        label="Size by"
        option={config.plotOptions.sizeBy}
        options={availableProperties}
        setOption={(sizeBy: string) => {
          setConfig({ ...config, plotOptions: { ...config.plotOptions, sizeBy } });
        }}
        availableNearestNeighbors={availableNearestNeighbors}
        availableClusters={availableClusters}
        availableProperties={availableProperties}
        extensions={[PlotSelectJaccardExtension, PlotSelectNNDiffExtension, PlotSelectEvalExtension]}
      />
      <PlotSelect
        inline={inline}
        id={`connectBy${id}Select`}
        label="Connect by"
        multi={true}
        option={config.plotOptions.connectBy}
        options={availableConnectByProperties}
        setOption={(option: string[]) => {
          setConfig({
            ...config,
            plotOptions: { ...config.plotOptions, connectBy: option.length === 0 ? null : option },
          });
        }}
        availableNearestNeighbors={availableNearestNeighbors}
        availableClusters={availableClusters}
        availableProperties={availableProperties}
        extensions={[PlotSelectJaccardExtension, PlotSelectNNDiffExtension, PlotSelectKNNExtension, PlotSelectClusterExtension, PlotSelectEvalExtension]}
      />
    </form>
  );
};
